#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 17 00:32:07 2022

@author: mtcontre
"""
import numpy as np
import pandas as pd
import geopandas
from itertools import compress

import os
import subprocess
from osgeo import osr
#from osgeo import gdal

def GetExtent(ds):
    """ Return list of corner coordinates from a gdal Dataset """
    xmin, xpixel, _, ymax, _, ypixel = ds.GetGeoTransform()
    width, height = ds.RasterXSize, ds.RasterYSize
    xmax = xmin + width * xpixel
    ymin = ymax + height * ypixel

    return (xmin, ymax), (xmax, ymax), (xmax, ymin), (xmin, ymin)

def ReprojectCoords(coords,src_srs,tgt_srs):
    """ Reproject a list of x,y coordinates. """
    trans_coords=[]
    transform = osr.CoordinateTransformation( src_srs, tgt_srs)
    for x,y in coords:
        x,y,z = transform.TransformPoint(x,y)
        trans_coords.append([x,y])
    return trans_coords


def convert_netcdf_to_geotiff(infile, outfile, band_number=1, verbose=False):
    # ds = gdal.Open(infile)
    # ds.SetProjection('WGS84')
    # met = ds.GetMetadata()
    # ulx = 
    # uly = 
    # lrx = 
    # lry = 
    with open(os.devnull, 'w') as devnull:
        return subprocess.check_call([
            "gdal_translate", infile, outfile,
            "-of", "GTiff",
            "-ot", "Float32",
            "-a_srs", "+proj=longlat +ellps=WGS84",
            #"-a_ullr", "-180","90","180","-90",
            "-a_ullr","-110","55","-89.95","39.99",
            "-b", str(band_number)

        ],
            stdout=None if verbose else devnull
        )


def readMesh(mesh):
    print('Reading mesh: '+ mesh)
    # Read nodes and connectivity table
    p = []
    t = []
    b = []
    with open(mesh + '.14','r') as msh:
        msh.readline()  #skip first line - comments
        aux = msh.readline() #skip first line - comments
        aux2 = aux.split()
        lp = int(aux2[1])
        lt = int(aux2[0])
        for count in range(lp+lt):
            if count<lp:
                items = msh.readline().split() 
                p.append([float(items[1]),float(items[2])])
                b.append(float(items[3]))
                count +=1
            elif count<lp+lt:
                items = msh.readline().split()
                t.append([int(items[2]),int(items[3]),int(items[4])])
                count +=1
        # Read boundary conditions
        
        # Open ocean boundary conditions
        items = msh.readline().split()
        nope = int(items[0]) # Number of open boundary conditions
        items = msh.readline().split()
        #neta = int(items[0]) # Number of open boundary nodes
        
        op = pd.DataFrame(data=[],columns = ['nvdll','nbdv'],\
                          dtype=object,index = range(nope))
        
        for nopes in range(nope):
            items = msh.readline().split()
            op.nvdll[nopes]=int(items[0]) # Number of nodes per open boundary condition
            nbdv = []
            for nbdvs in range(op.nvdll[nopes]):
                items = msh.readline().split()
                nbdv.append(int(items[0]))
            op.nbdv[nopes] = nbdv
        
        # Land boundary conditions
        items = msh.readline().split()
        nbou = int(items[0]) # Number of land boundary conditions
        items = msh.readline().split()
        #nvel = int(items[0]) # Number of land boundary nodes
        
        bd = pd.DataFrame(data=[],\
                          columns = ['ibtype','nvell','barincfsb','barincfsp','nbvv','ibconn','barinht'],\
                          dtype=object,index = range(nbou))
        for nbous in range(nbou):
            items = msh.readline().split()
            bd.nvell[nbous] = int(items[0])
            bd.ibtype[nbous] = int(items[1])
            
            if bd.ibtype[nbous]==24: # Levees boundary conditions
                nbvv = []; ibconn = []; barinht = []; barincfsb = []; barincfsp = [];
                for nvells in range(bd.nvell[nbous]):
                    items = msh.readline().split()
                    nbvv.append(int(items[0]))
                    ibconn.append(int(items[1]))
                    barinht.append(float(items[2]))
                    barincfsb.append(int(float(items[3])))
                    barincfsp.append(int(float(items[4])))
                bd.nbvv[nbous] = nbvv
                bd.ibconn[nbous] = ibconn
                bd.barinht[nbous] = barinht
                bd.barincfsb[nbous] = barincfsb
                bd.barincfsp[nbous] = barincfsp
            
            if bd.ibtype[nbous]== 20 or bd.ibtype[nbous]==21 or bd.ibtype[nbous]==22: # Land boundary conditions
                nbvv = []; 
                for nvells in range(bd.nvell[nbous]):
                    items = msh.readline().split()
                    nbvv.append(int(items[0]))
                bd.nbvv[nbous] = nbvv

    msh.close()

    # Transform lists to arrays
    p = np.array(p)
    t = np.array(t)
    #b = np.array(b)

    from sklearn.neighbors import KDTree
    tree_p = KDTree(p)


    from shapely.geometry import Point

    geoms = [Point(p[i]) for i in range(len(p))] #geoms[0].coords.xy
    data_dict = {'b':b}
    p = geopandas.GeoDataFrame(data_dict, geometry=geoms, crs="EPSG:4326")
    
    print('Mesh succesfully read')
    
    return p,t,tree_p,bd

def PropSecNat(coords,Z):
    
    from scipy.spatial import distance

    # Funcion que calcula las propiedades geometricas de cualquier seccion
    # [A Pm B Rh Dh]=PropSecNat(seccion, h)

    # Inputs:
    # "seccion" debe ser una matriz de coordenadas de dos columnas
    # que contiene los puntos que definen la matriz
    # La unica condicion es que los puntos se den en el orden de recorrido
    # de la forma de la seccion
    # "h" es la altura de agua sobre el punto mas bajo de la seccion

    # Outputs:
    # A = Área
    # B = Ancho superficial
    # Pm = Perímetro mojado
    # Rh = Radio hidráulico (A/Pm)
    # Dh = Profundidad hidráulica(A/B)
    
    
    aux = distance.pdist(coords,'euclidean')
    aux = np.radians(aux[0:len(coords)-1])*6371*1000
    Y = np.hstack([0,aux])
    # if len(aux)>1:
    #     Y = np.hstack([0,aux[0],aux+aux[0],aux[-1]+aux[0]+(aux[-1]-aux[-2])])
    # else:
    #     Y = np.hstack([0,aux[0],aux+aux[0],aux[-1]+aux[0]+(aux[-1])])
    
    n=len(Y)
    Z = -1*Z
    #Z = np.hstack([0.3, Z, 0.3])

    #Z0=min(Z)
    Zh=0
    H=-1*Z
    #j=find(Z>Zh)
    #j = Z[Z>Zh]


    # I: funcion indicatriz
    I=np.full(n,True)
    I[Z>Zh]=False


    A = np.zeros(n-1);Pm =np.zeros(n-1);B=np.zeros(n-1);CG=np.zeros(n-1)
    #vamos por cada par de nodos
    for i in range(n-1):
        if ~I[i]:    #primer nodo seco
            if ~I[i+1]:  #y segundo nodo seco (ambos fuera)
                A[i]=0
                Pm[i]=0
                B[i]=0
                CG[i]=0
            else:        #y segundo nodo mojado (CC:BB izq)
                Yizq=Y[i]-(Y[i+1]-Y[i])/(H[i+1]-H[i])*H[i]
                B[i]=Y[i+1]-Yizq
                A[i]=H[i+1]*B[i]/2
                Pm[i]=np.sqrt(H[i+1]**2+B[i]**2)
                CG[i]=A[i]*H[i+1]/3
        else:        #primer nodo mojado
           if ~I[i+1]:  #y segundo nodo seco (CC:BB der)
                Yder=Y[i]-(Y[i+1]-Y[i])/(H[i+1]-H[i])*H[i]
                B[i]=Yder-Y[i]
                A[i]=H[i]*B[i]/2
                Pm[i]=np.sqrt(H[i]**2+B[i]**2)
                CG[i]=A[i]*H[i]/3
           else:        #y segundo nodo mojado (nodo central)
                B[i]=Y[i+1]-Y[i]
                A[i]=(H[i]+H[i+1])/2*B[i]
                Pm[i]=np.sqrt(B[i]**2+(H[i]-H[i+1])**2)
                CG[i]=A[i]*(H[i]+H[i+1])/2/2
  
  
  
    #A=sum(A)
    #Pm=sum(Pm)
    #if I[0]==1: Pm=Pm+H[0]

    #if I[n-1]==1: Pm=Pm+H[n-1]
   
    #Rh=A/Pm
    #B=sum(B)
    #Dh=A/B
    #CG=sum(CG)/A
    
    
    
    return A,B


def getFeatures(gdf):
    """Function to parse features from GeoDataFrame in such a manner that rasterio wants them"""
    import json
    return [json.loads(gdf.to_json())['features'][0]['geometry']]

def fixmesh(p,t,ptol=1024*np.finfo(float).eps):
#FIXMESH  Remove duplicated/unused nodes and fix element orientation.
#   [P,T]=FIXMESH(P,T)

#   Copyright (C) 2004-2012 Per-Olof Persson. See COPYRIGHT.TXT for details.
    

    snap = np.amax(np.amax(p, axis=0)-np.amin(p, axis=0),axis=0)*ptol
    aux,ix,jx = np.unique(np.round(p/snap)*snap,return_index=True,return_inverse=True,axis = 0)
    # jx = sorted(ix)
    #indices = np.argsort(ix)
    ix_sorted = np.sort(ix)
    
    numRows = len(p[:,0])
    indSortA = np.argsort(p[:,0],axis=0)
    sortA= p[indSortA,:]
    
    groupsSortA = sortA[0:numRows-1,:] != sortA[1:numRows,:]
    groupsSortA = np.any(groupsSortA,axis=1)
    groupsSortA = np.concatenate((np.array([True]),groupsSortA))
    invIndSortA = indSortA
    invIndSortA[invIndSortA] = np.linspace(0,numRows-1,numRows)  # Find inverse permutation.   
    logIndA = groupsSortA[invIndSortA]   # Create new logical by indexing into groupsSortA.  
    c = p[logIndA,:]
        
    indSortC = np.argsort(c[:,0],axis=0)    # Sort C to get index.                  
    lengthGroupsSortA = np.diff(list(compress(range(len((groupsSortA, True))), (groupsSortA, True)))) # Determine how many of each of the above indices there are in IC.              
    diffIndSortC = np.diff(indSortC)
    diffIndSortC = (indSortC[0], diffIndSortC)
                    
    indLengthGroupsSortA = np.cumsum((1, lengthGroupsSortA));  # Get the correct amount of each index.
    indLengthGroupsSortA[-1] = []
                    
    indCOrderedBySortA = np.zeros((numRows,0))
    indCOrderedBySortA[indLengthGroupsSortA,0] = diffIndSortC  # Since indCOrderedBySortA is not already established as a column,                    
    indCOrderedBySortA = np.cumsum(indCOrderedBySortA)
    jx = indCOrderedBySortA[invIndSortA]               

       
    p = p[ix,:]
    t = np.reshape(jx[t],(np.size(t,0),np.size(t,1)))
    
    pix,jx1=np.unique(t,return_inverse=True,return_index=False)
    t = np.reshape(jx1,(np.size(t,0),np.size(t,1)))
    p = p[pix,:]
    pix=ix[pix]
    
    if np.size(t,1) == np.size(p,1)+1:
        flipped = (OM.simpvol(p,t)<0)
        aux = t[flipped,1]
        aux2 = t[flipped,0]
        t[flipped,0] = aux
        t[flipped,1] = aux2
    


    return p,t,pix


def simpvol(p,t):
#SIMPVOL Simplex volume.
#  V=SIMPVOL(P,T)

#  Copyright (C) 2004-2012 Per-Olof Persson. See COPYRIGHT.TXT for details.
    import numpy as np
    aux = np.size(p,1)
    if aux ==1:
      d12 = p[t[:,1],:]-p[t[:,0],:]
      v = d12;
    elif aux ==2:
      d12 = p[t[:,1],:]-p[t[:,0],:]
      d13 = p[t[:,2],:]-p[t[:,0],:]
      v = (d12[:,0]*d13[:,1]-d12[:,1]*d13[:,0])/2
    elif aux ==3:
      d12 = p[t[:,1],:]-p[t[:,0],:]
      d13 = p[t[:,2],:]-p[t[:,0],:]
      d14 = p[t[:,3],:]-p[t[:,0],:]
      v = np.dot(np.cross(d12,d13,2),d14,2)/6
    else:
      v = np.zeros(np.size(t,0),1)
      for ii in range(np.size(t,0)):
        A = np.zeros(np.size(p,1)+1)
        A[:,0] = 1
        for jj in range(np.size(p,1)+1):
          A[jj,1:-1] = p[t[ii,jj],:]
        v[ii] = np.linalg.det(A)
      v = v/np.math.factorial(np.size(p,1))

    return v




def renum(p,t,b=np.empty(0),bx=np.empty(0),by=np.empty(0)):
            
    import numpy as np
    import scipy.sparse  as sp
    from scipy.sparse.csgraph import reverse_cuthill_mckee


    np1=len(p[:,0])
    nt=len(t[:,0])
    #Calculate adjacency matrix of t
    aux0 = np.ones((6*nt,1)).flatten(); aux0 = aux0.astype(np.int64);
    aux1 = t[:,(0,0,1,1,2,2)].flatten()-1
    aux2 = t[:,(1,2,0,2,0,1)].flatten()-1
    S = sp.csr_matrix((aux0,(aux1,aux2)),shape=(np1,np1))
  
    W = np.sum(S,axis=1)
    try:
        any(W==0)
    except Exception:
        print("ERROR: Invalid mesh. Hanging nodes found. Retriangulate.")
        
            
    # calc bw
    i,j,v = sp.find(S)
    bw = max(i-j) + 1
    print('Initial bandwidth is ',str(bw))
            
    #renumber with minimizing bw
    perm = reverse_cuthill_mckee(S,symmetric_mode=False)
    perm = eng.symrcm((S),nargout=1)
    R    = S[perm]
    prn  = p[perm,:]
    if b.size>0:
        brn  = b[perm]
        b    = brn
        
    if bx.size>0:
        brn  = bx[perm]
        bx    = brn
    if by.size>0:            
        brn  = by[perm,:]
        by    = brn
            
    p = prn
    perm_inv = np.zeros(np1)
    perm_inv[perm[0:np1]] = np.linspace(1,np1,num=np1)
    ttemp = t 
    for ie in range(nt):
        nm1 = ttemp[ie,0]
        nm2 = ttemp[ie,1]
        nm3 = ttemp[ie,2]
        ttemp[ie,0] = perm_inv[nm1]
        ttemp[ie,1] = perm_inv[nm2]
        ttemp[ie,2] = perm_inv[nm3]
            
    t = ttemp 
    # compare bw
    i,j,v = sp.find(R)
    bw = max(i-j) + 1
    print("Renumbered bandwidth is ",str(bw))
            
           
    return p,t


# def obj = Make_Mesh_Boundaries_Traversable( obj, dj_cutoff, nscreen, proj )
# %  obj =  Make_Mesh_Boundaries_Traversable(obj,dj_cutoff,nscreen)
# %  A msh object (containing p and t) is  "cleaned" and returned.
# %  ncscreen ~= 0 will display info to screen
# %  See "Exterior Check" description below for definition of dj_cutoff
# %
# %  Alternates between checking interior and exterior portions
# %  of the graph exhaustively until convergence is obtained, defined as:
# %  Having no nodes connected to more than 2 boundary edges.
# %
# %  Interior Check: Deletes elements that are within the interior of the
# %  mesh so that no nodes are connected to more than 2 boundary edges. For
# %  example, a spit could become very thin in a middle portion so that you
# %  a node is connected to two elements but four boundary edges, in a
# %  bow-tie type formation. This code will delete one of those connecting
# %  elements to ensure the spit is continous and only two boundary edges
# %  are connected to that node. In the case of a choice between elements to
# %  delete, the one with the lowest quality is chosen.
# %
# %  Exterior Check: Finds small disjoint portions of the graph and removes
# %  them using a breadth-first search. The individual disjoint portions are
# %  removed based on dj_cutoff.
# %  dj_cutoff <= 1 indicates the proportional of the total area of the mesh
# %  that, below which, a disjoint portion is discarded.
# %  dj_cutoff > 1 indicates the absolute area in km2 that, below which, a
# %  disjoint portion is discarded.
# %  dj_cutoff can be used to keep only the largest area of the mesh
# %  (e.g. dj_cutoff = 0.5) or a few largest areas (e.g. dj_cutoff = 0.1),
# %  or it can be set to ensure that lakes or polders of certain sizes are
# %  kept (e.g. dj_cutoff = 1 [km2])
# %
# %  Copyright (C) 2018  Keith Roberts & William Pringle
# %  Algorithm written by William Pringle, CHL, UND 2017
# %  Improvements by Keith Roberts, CHL, UND 2017
# %  WJP: Organisation of code Jan 2018
# %  KJR: Speed-up, CHL, UND 2018.
# %
# %  This program is free software: you can redistribute it and/or modify
# %  it under the terms of the GNU General Public License as published by
# %  the Free Software Foundation, either version 3 of the License, or
# %  (at your option) any later version.
# %
# %  This program is distributed in the hope that it will be useful,
# %  but WITHOUT ANY WARRANTY; without even the implied warranty of
# %  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# %  GNU General Public License for more details.
# %
# %  You should have received a copy of the GNU General Public License
# %  along with this program.  If not, see <http://www.gnu.org/licenses/>.

# % Entering the code
# disp('Making mesh boundaries traversable...');
# % This is to avoid heaps of warnings in the triangulation call which are
# % unneccesary
# warning('off','all')

# if nargin < 3
#     % display outputs
#     nscreen = 1;
# end
# if nargin < 4
#     % already projected
#     proj = 1;
# end

# % Get p and t out of obj
# p = obj.p; t = obj.t;

# % Delete disjoint nodes
# [p,t,xf] = fixmesh(p,t);
# if ~isempty(obj.b)
#     obj.b = obj.b(xf);
# end
# if ~isempty(obj.bx)
#     obj.bx = obj.bx(xf);
# end
# if ~isempty(obj.by)
#     obj.by = obj.by(xf);
# end

# % Get boundary edges and nodes
# [etbv,vxe] = extdom_edges2( t, p ) ;

# % WJP sometimes we may wanna delete some exterior portions even with a valid
# % mesh so allow entry even in this case.
# if dj_cutoff > 0 && numel(etbv) == numel(vxe)
#     etbv(end+1,:) = 1;
# end
# % Loop until all the nodes only have two boundary edges
# %(the number of boundary edges will equal number of boundary nodes)

# while numel(etbv) > numel(vxe)
    
#     % Delete elements in the exterior of the mesh
#     t = delete_exterior_elements(p,t,dj_cutoff,nscreen,proj);
    
#     % Delete disjoint nodes
#     [p,t,xf] = fixmesh(p,t);
#     if ~isempty(obj.b)
#         obj.b = obj.b(xf);
#     end
#     if ~isempty(obj.bx)
#         obj.bx = obj.bx(xf);
#     end
#     if ~isempty(obj.by)
#         obj.by = obj.by(xf);
#     end
    
#     % Delete elements in the interior of the mesh
#     t = delete_interior_elements(p,t,nscreen);
    
#     % Delete disjoint nodes
#     [p,t,xf] = fixmesh(p,t);
#     if ~isempty(obj.b)
#         obj.b = obj.b(xf);
#     end
#     if ~isempty(obj.bx)
#         obj.bx = obj.bx(xf);
#     end
#     if ~isempty(obj.by)
#         obj.by = obj.by(xf);
#     end
    
#     % Get boundary edges and nodes
#     [etbv,vxe] = extdom_edges2( t, p ) ;
    
#     %if numel(vxe) > numel(etbv)
#     %   error(['number of boundary vertices larger than boundary edges', ...
#     %          ', try a larger dj_cutoff'])
#     %end
# end
# % Finished cleaning
# disp('ALERT: finished cleaning up mesh..');
# % Turn warnings back on
# warning('on','all')
# % Put back into the msh object
# obj.p = p; obj.t = t;
# end
# % The sub-functions...
# %% Delete elements outside the main mesh depending on dj_cutoff input
# % dj_cutoff >= 1
# %    area in km2
# % dj_cutoff < 1
# %    proportion of the total mesh area
# function t = delete_exterior_elements(p,t,dj_cutoff,nscreen,proj)
# %
# if dj_cutoff == 0
#     if nscreen
#         disp('dj_cutoff is zero; do nothing in delete_exterior_elements')
#     end
#     % do nothing
#     return;
# end
# L = size(t,1);
# t1 = t; t = [];
# if proj
#     % has already been projected so need to convert back to lat-lon
#     [X,Y] = m_xy2ll(p(:,1),p(:,2));
# else
#     % not projected so keep the lat-lon;
#     X = p(:,1); Y = p(:,2);
# end
# % calculate area
# A = sum(polyarea(X(t1)',Y(t1)').*cosd(mean(Y(t1)')));
# %A = sum(polyarea(X(t1(:,1:3))',Y(t1(:,1:3))'));
# An = A;
# if dj_cutoff >= 1
#     Re2 = 111^2; % Re2 = (6378.137)^2;
#     An = Re2*An;
#     % Absolute area
#     while An > dj_cutoff
        
#         % Peform the Breadth-First-Search to get nflag
#         nflag = BFS(p,t1);
        
#         % Get new triangulation and its area
#         t2 = t1(nflag == 1,:);
#         An = Re2*sum(polyarea(X(t2)',Y(t2)').*cosd(mean(Y(t2)')));
        
#         % If large enough at t2 to the triangulation
#         if An > dj_cutoff
#             t = [t; t2];
#         end
#         % Delete where nflag == 1 since this patch didn't meet the fraction
#         % limit criterion.
#         t1(nflag == 1,:) = [];
#         % Calculate the remaining area
#         An = Re2*sum(polyarea(X(t1)',Y(t1)').*cosd(mean(Y(t1)')));
        
#     end
# elseif dj_cutoff > 0
#     % Area proportion
#     while An/A > dj_cutoff
        
#         % Peform the Breadth-First-Search to get nflag
#         nflag = BFS(p,t1);
        
#         % Get new triangulation and its area
#         t2 = t1(nflag == 1,:);
#         An = sum(polyarea(X(t2)',Y(t2)').*cosd(mean(Y(t2)')));
        
#         % If large enough at t2 to the triangulation
#         if An/A > dj_cutoff
#             t = [t; t2];
#         end
#         % Delete where nflag == 1 since this patch didn't meet the fraction
#         % limit criterion.
#         t1(nflag == 1,:) = [];
#         % Calculate the remaining area
#         An = sum(polyarea(X(t1)',Y(t1)').*cosd(mean(Y(t1)')));
        
#     end
# elseif dj_cutoff < 0
#     %error('Keep cannot be negative 0')
#     tic
#     while length(t1(:,1))>=-1*dj_cutoff
        
#         % Peform the Breadth-First-Search to get nflag
#         nflag = BFS(p,t1);
        
#         % Get new triangulation and its area
#         t2 = t1(nflag == 1,:);
        
#         % If large enough at t2 to the triangulation
#         if length(t2(:,1)) > -1*dj_cutoff
#             t = [t; t2];
#         end
#         % Delete where nflag == 1 since this patch didn't meet the fraction
#         % limit criterion.
#         t1(nflag == 1,:) = [];
#         [length(t2(:,1)),length(t1(:,1))]
        
#     end
#     toc
# else
#     % keep is zero, do nothing
# end
# if nscreen
#     disp(['  ACCEPTED: deleting ' num2str(L-size(t,1)) ...
#         ' elements outside main mesh']) ;
# end
# if size(t,1) < 1
#     error(['All elements have been deleted... something wrong? ' ...
#         'dj_cutoff is set to' num2str(dj_cutoff)])
# end

# end

# %% Delete interior elements
# function t = delete_interior_elements(p,t,nscreen)
# % Get updated boundary edges.
# etbv = extdom_edges2( t, p ) ;
# % Get all nodes that are on edges.
# [nodes_on_edge,~,n] = unique(etbv(:));
# % Count how many edges a node appears in.
# I     = accumarray(n,1:numel(n),[],@(x){x});
# count = cellfun('length',I);
# %
# [vtoe,nne] = VertToEle(t);
# % Get the nodes which appear more than twice and delete element connected
# % to these nodes where all nodes of element are on boundary edges
# del_elem_idx = [];
# for i = nodes_on_edge(count > 2)'
#     con_elem = vtoe(1:nne(i),i);
#     n = 0; del_elem = [];
#     for elem = con_elem'
#         I = etbv(:) == t(elem,1);
#         J = etbv(:) == t(elem,2);
#         K = etbv(:) == t(elem,3);
#         % All nodes on element are boundary edges
#         if any(I) && any(J) && any(K)
#             n = n + 1;
#             del_elem(n) = elem;
#         end
#     end
#     if n == 1
#         % Only one element to delete.
#         del_elem_idx(end+1) = del_elem;
#     elseif n > 1
#         % Delete worst quality qualifying element.
#         tq = gettrimeshquan( p, t(del_elem,:));
#         [~,idx] = min(tq.qm);
#         del_elem_idx(end+1) = del_elem(idx);
#     else
#         % No connected elements have all nodes on boundary edge so we
#         % select the worst quality connecting element.
#         tq = gettrimeshquan( p, t(con_elem,:));
#         [~,idx] = min(tq.qm);
#         del_elem_idx(end+1) = con_elem(idx);
#     end
# end

# if nscreen
#     disp(['  ACCEPTED: deleting ' num2str(length(del_elem_idx)) ...
#         ' elements inside main mesh'])
# end
# t(del_elem_idx,:) = [];

# end

# function nflag = BFS(p,t1)

# % Select a random element.
# EToS =  randi(size(t1,1),1);

# % Get element-to-element connectivity.
# tri = triangulation(t1,p);
# nei = tri.neighbors;

# % Traverse grid deleting elements outside.
# ic = zeros(ceil(sqrt(size(t1,1))*2),1);
# ic0 = zeros(ceil(sqrt(size(t1,1))*2),1);
# nflag = zeros(size(t1,1),1);
# ic(1) = EToS;
# icc  = 1;

# % Using BFS loop over until convergence is reached (i.e., we
# % obtained a connected region).
# while icc
#     ic0(1:icc) = ic(1:icc);
#     icc0 = icc;
#     icc = 0;
#     for nn = 1:icc0
#         i = ic0(nn);
#         % Flag the current element as OK
#         nflag(i) = 1;
#         % Search neighbouring elements
#         nb = nei(i,:);
#         nb(isnan(nb)) = [];
#         % Flag connected neighbors as OK
#         for nnb = 1:length(nb)
#             if ~nflag(nb(nnb))
#                 icc = icc + 1;
#                 ic(icc) = nb(nnb);
#                 nflag(nb(nnb)) = 1;
#             end
#         end
#     end
# end
# end
