# laminatelib.py
import numpy as np


# example_layup = [ {'mat':m1 , 'ori':  0  , 'thi':0.5} ,
#            {'mat':m1 , 'ori': 90  , 'thi':0.5} ,
#            {'mat':m1 , 'ori': 45  , 'thi':0.5} ,
#            {'mat':m1 , 'ori':-45  , 'thi':0.5}  ]
#





def S2D(m):
    return np.array([[        1/m['E1'], -m['v12']/m['E1'],          0],
                     [-m['v12']/m['E1'],         1/m['E2'],          0],
                     [                0,                 0, 1/m['G12']]], float)

def Q2D(m):
    S=S2D(m)
    return np.linalg.inv(S)



def T2Ds(a):
    c,s = np.cos(np.radians(a)), np.sin(np.radians(a))
    return np.array([[ c*c ,  s*s ,   2*c*s],
                     [ s*s ,  c*c ,  -2*c*s],
                     [-c*s,   c*s , c*c-s*s]], float)

def T2De(a):
    c,s = np.cos(np.radians(a)), np.sin(np.radians(a))
    return np.array([[   c*c,   s*s,     c*s ],
                     [   s*s,   c*c,    -c*s ],
                     [-2*c*s, 2*c*s, c*c-s*s ]], float)

def Q2Dtransform(Q,a):
    return np.dot(np.linalg.inv( T2Ds(a) ) , np.dot(Q,T2De(a)) )

def laminateThickness(layup):
    return sum([layer['thi'] for layer in layup])



def laminateStiffnessMatrix(layup):
    ABD=np.zeros((6,6),float)
    hbot = -laminateThickness(layup)/2        # bottom of first layer
    for layer in layup:
        Qt = Q2Dtransform(Q2D(layer['mat']), layer['ori'])
        htop = hbot + layer['thi']   # top of current layer
        ABD[0:3,0:3] = ABD[0:3,0:3] +       Qt*(htop   -hbot   )
        ABD[0:3,3:6] = ABD[0:3,3:6] + (1/2)*Qt*(htop**2-hbot**2)
        ABD[3:6,0:3] = ABD[3:6,0:3] + (1/2)*Qt*(htop**2-hbot**2)
        ABD[3:6,3:6] = ABD[3:6,3:6] + (1/3)*Qt*(htop**3-hbot**3)
        hbot=htop                    # for the next layer
    return ABD

def solveLaminateLoadCase(ABD,**kwargs):
    # ABD: the laminate stiffness matrix
    # kwargs: prescribed loads or deformations
    loads=np.zeros((6))
    defor=np.zeros((6))
    loadKeys=('Nx','Ny','Nxy','Mx','My','Mxy')   # valid load keys
    defoKeys=('ex','ey','exy','Kx','Ky','Kxy')   # valid deformation keys
    knownKeys=['Nx','Ny','Nxy','Mx','My','Mxy']  # initially assume known loads
    knownVals=np.array([0,0,0,0,0,0],float)
    for key in kwargs:
        if key in loadKeys:
            i=loadKeys.index(key)
            knownKeys[i]=key
            knownVals[i]=kwargs[key]
        if key in defoKeys:
            i=defoKeys.index(key)
            knownKeys[i]=key
            knownVals[i]=kwargs[key]
    M1=-np.copy(ABD)
    M2= np.copy(ABD)
    ID= np.identity(6)
    for k in range(0,6):
        if knownKeys[k] in loadKeys:
            M2[:,k]=-ID[:,k]
        else:
            M1[:,k]=ID[:,k]
    sol=np.dot(  np.linalg.inv(M1),   np.dot(M2,knownVals))
    for i in range(0,6):
        if knownKeys[i] == loadKeys[i]:
            loads[i]=knownVals[i]
            defor[i]=sol[i]
        if knownKeys[i] == defoKeys[i]:
            defor[i]=knownVals[i]
            loads[i]=sol[i]
    return loads,defor

def thermalLoad(layup,dT):
    NTMT=np.zeros(6,float) #thermal loads
    hbot = -laminateThickness(layup)/2
    for layer in layup:
        m = layer['mat']
        Q = Q2D(m)
        a = layer['ori']
        Qt = Q2Dtransform(Q,a)
        a123=[m['a1'], m['a2'], 0]
        aXYZ=np.dot( T2De(-layer['ori']) , a123 )
        htop = hbot + layer['thi']
        NTMT[0:3] = NTMT[0:3] +     dT*(  np.dot(Qt,aXYZ)    )*(htop   -hbot   )
        NTMT[3:6] = NTMT[3:6] + 0.5*dT*(  np.dot(Qt,aXYZ)    )*(htop**2-hbot**2)
        hbot=htop
    return NTMT

def layerResults(layup,deformations):
    # An empty list that shall be populated with results from all layers:
    results=[]
    # The bottom coordinate of the laminate:
    bot = -laminateThickness(layup)/2
    for layer in layup:
        # the top of current layer is just the bottom + the thickness:
        top=bot+layer['thi']
        # creating a dictionary with the two coordinates:
        h={'bot':bot, 'top':top}
        # computing the strains in laminate coordinate system:
        strnXYZbot=deformations[0:3]+bot*deformations[3:6]
        strnXYZtop=deformations[0:3]+top*deformations[3:6]
        # put the strains into a dictionary:
        strnXYZ={'bot':strnXYZbot, 'top':strnXYZtop}
        # the transformation matrix for strains:
        Te=T2De(layer['ori'])
        # transforming the strains into the material coordinate system:
        strn123bot=np.dot(Te, strnXYZbot)
        strn123top=np.dot(Te, strnXYZtop)
        # the strains at top and bottom as a dictionary:
        strn123={'bot':strn123bot, 'top':strn123top}
        # all strains as a new dictionary:
        strn={'xyz':strnXYZ, '123':strn123}
        # stiffness matrix of the material:
        Q=Q2D(layer['mat'])
        # Hooke's law: finding the stresses in material system:
        strs123bot=np.dot(Q,strn123bot)
        strs123top=np.dot(Q,strn123top)
        # the material stresses as a dictionary having bottom and topp results:
        strs123={'bot':strs123bot, 'top':strs123top}
        # transformation matrix for stress, negativ rotation applies now,
        # since this is from material system to laminat system (reversed..)
        Ts=T2Ds(-layer['ori'])
        # transforming stresses into laminate system:
        strsXYZbot=np.dot(Ts,strs123bot)
        strsXYZtop=np.dot(Ts,strs123top)
        # organizing the stresses for the two location in a dictionary
        strsXYZ={'bot':strsXYZbot, 'top':strsXYZtop}
        # all stresses as a new nested dictionary:
        strs={'xyz':strsXYZ, '123':strs123}
        # Failure criteria...follows the same logic as above
        MSbot=fE2DMS(strs123bot,layer['mat'])
        MStop=fE2DMS(strs123top,layer['mat'])
        MS={'bot':MSbot, 'top':MStop}
        MEbot=fE2DME(strs123bot,layer['mat'])
        MEtop=fE2DME(strs123top,layer['mat'])
        ME={'bot':MEbot, 'top':MEtop}
        TWbot=fE2DTW(strs123bot,layer['mat'])
        TWtop=fE2DTW(strs123top,layer['mat'])
        TW={'bot':TWbot, 'top':TWtop}
        fail={'MS':MS, 'ME':ME, 'TW':TW}
        # and finally put everything into a top level dictionary,
        # and add that to the list
        results.append( {'h':h , 'strain':strn, 'stress':strs, 'fail':fail }  )
        bot=top
    return results