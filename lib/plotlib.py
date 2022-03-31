#plotlib.py

import numpy as np
import matplotlib.pyplot as plt

def plotMatIndex(names,xdata,ydata,xlabel,ylabel):
    import numpy as np
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(10,6))
    area = 50
    for k in range(0,len(names)):
        ax.text(xdata[k]*1.01,ydata[k]*1.01,names[k])
    ax.grid(True)
    ax.scatter(xdata, ydata, s=area, c='green', alpha=0.5)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    plt.show()


    ##EXAMPLE:
    #names=['E-glass', 'T1100', 'M60', 'Steel','Cryptonite']
    #xdata=[2.55, 1.79, 1.93, 7.8, 5.0]
    #ydata=[76, 324, 588, 201,400]
    #%matplotlib inline
    #plotMatIndex(names,xdata,ydata,'Density','Modulus')


def plotSquareArrayOfFibers(Vf):
    import matplotlib.pyplot as plt
    import numpy as np
    fig, ax = plt.subplots(figsize=(3, 3))
    a = 1.0
    ax.set_xlim(-2.1 * a, 2.1 * a)
    ax.set_ylim(-2.1 * a, 2.1 * a)
    ax.set_axis_off()

    d = ((4 * Vf * a ** 2) / (np.pi)) ** 0.5
    circle1 = plt.Circle((-a, -a), d, color='black', fc='silver')
    circle2 = plt.Circle((a, -a), d, color='black', fc='silver')
    circle3 = plt.Circle((a, a), d, color='black', fc='silver')
    circle4 = plt.Circle((-a, a), d, color='black', fc='silver')
    ax.add_artist(circle1)
    ax.add_artist(circle2)
    ax.add_artist(circle3)
    ax.add_artist(circle4)
    ax.plot((-a, a, a, -a, -a), (-a, -a, a, a, -a), '--', color='black',
            linewidth=1)
    plt.show()

    # EXAMPLE:
    #% matplotlib
    #inline
    #plotSquareArrayOfFibers(Vf=0.55)


def plotRandomDistributionOfFibers(dx, dy, d, n):
    # dx and dy are dimensions
    # d is a list or 1D array of diameters
    # n is the number of attempts
    import numpy as np
    rx = dx * np.random.rand(n)
    ry = dy * np.random.rand(n)
    x, y = [], []  # lists of coordinates verified to fit
    k = 0  # counter: no. of fibers already added

    for j in range(0, n):  # all randomly generated coordinates
        fits = True  # initially assuming that the fiber fits
        for p in range(0, k):  # all previously verified coordinates
            distance = ((rx[j] - x[p]) ** 2 + (
                        ry[j] - y[p]) ** 2) ** 0.5  # vector length
            if distance < (d[k] + d[
                p]) / 2:  # does not fit if the distance is less than average diameter of the two
                fits = False
        if fits == True:  # append only the ones that fit
            x.append(rx[j])
            y.append(ry[j])
            k = k + 1
    print('Number of fibers:', k)
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.set_xlim(0, dx)
    ax.set_ylim(0, dy)
    ax.set_axis_off()
    for i in range(0, len(x)):
        circle1 = plt.Circle((x[i], y[i]), d[i] / 2, color='black', fc='silver')
        ax.add_artist(circle1)


    ## EXAMPLE
    #import numpy as np

    #
    #fibd = np.random.normal(10, 1, 1000)  # (mean, std, count)
    #% matplotlib
    #inline
    #plotRandomDistributionOfFibers(dx=200, dy=200, d=fibd, n=20000)

def plotDistribution(d):
    import matplotlib.pyplot as plt
    fig,ax = plt.subplots(figsize=(8,6))
    n,bins,patches=plt.hist(fibd, 20, density=True,facecolor='teal', alpha=0.6 )
    ax.set_xlabel('Fiber diameter(um)')
    ax.set_ylabel('Normalized number of fibers')
    ax.set_title('Fiber diameter distribution')
    ax.set_xlim(0,)
    ax.grid(True)
    plt.show()

    ##EXAMPLE:
    #import numpy as np
    #fibd = np.random.normal(10, 1, 1000)
    #%matplotlib inline
    #plotDistribution(fibd)

def illustrateStrains(ex,ey,exy,scaleFactor):
    import numpy
    import matplotlib.pyplot as plt
    x=numpy.array([-0.5,0.5,0.5,-0.5,-0.5])
    y=numpy.array([-0.5,-0.5,0.5,0.5,-0.5])
    ux=(ex*x+0.5*exy*y)*scaleFactor
    uy=(ey*y+0.5*exy*x)*scaleFactor
    xd=x+ux
    yd=y+uy
    fig,ax=plt.subplots(figsize=(4,4))
    ax.set_xlim(-1,1)
    ax.set_ylim(-1,1)
    ax.set_axis_off()
    ax.set_title('Scale factor='+str(scaleFactor))
    text=r'$\varepsilon_x=$'+str(ex)+'\n'+r'$\varepsilon_y=$'+str(ey)+'\n'+r'$\gamma_{xy}=$'+str(exy)
    ax.text(0,0,text, ha='center',va='center')
    ax.plot(x,y,'--',color='lightgreen')
    ax.plot(xd,yd,'-',color='darkgreen')
    ax.arrow(-0.9, -0.9, 1.0, 0.0, head_width=0.05, head_length=0.05, fc='black', ec='black')
    ax.arrow(-0.9, -0.9, 0.0, 1.0, head_width=0.05, head_length=0.05, fc='black', ec='black')
    ax.text(0.2,-0.9,'x')
    ax.text(-0.9,0.2,'y')

    ## EXAMPLE:
    #%matplotlib inline
    #illustrateStrains( ex=0.003, ey=-0.001, exy=-0.003, scaleFactor=100)


def illustrateLayup(layup, size=(4, 4)):
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    fig, ax = plt.subplots(figsize=size)
    tot = 0
    for layer in layup:
        tot = tot + layer['thi']
    hb = -tot / 2.0
    for layer in layup:
        ht = hb + layer['thi']
        if layer['ori'] > 0:
            fco = 'lightskyblue'
        if layer['ori'] < 0:
            fco = 'pink'
        if layer['ori'] == 0:
            fco = 'linen'
        if layer['ori'] == 90:
            fco = 'silver'
        p = Rectangle((-0.6, hb), 1.2, layer['thi'], fill=True,
                      clip_on=False, ec='black', fc=fco)
        ax.add_patch(p)
        mid = (ht + hb) / 2.0
        ax.text(0.62, mid, str(layer['ori']), va='center')
        hb = ht
    ax.set_xlim(-1, 1)
    ax.set_ylim(-1.1 * tot / 2.0, 1.1 * tot / 2.0)
    ax.get_xaxis().set_visible(False)
    ax.plot((-1, -0.8), (0, 0), '--', color='black')
    ax.plot((0.8, 1.0), (0, 0), '--', color='black')
    plt.show()

    # # EXAMPLE:
    # % matplotlib
    # inline
    # m1 = {'E1': 40000, 'E2': 10000, 'v12': 0.3, 'G12': 3000}
    # layup1 = [{'mat': m1, 'ori': 0, 'thi': 0.5},
    #           {'mat': m1, 'ori': 90, 'thi': 0.5},
    #           {'mat': m1, 'ori': 45, 'thi': 0.5},
    #           {'mat': m1, 'ori': -45, 'thi': 0.5},
    #           {'mat': m1, 'ori': -45, 'thi': 0.5},
    #           {'mat': m1, 'ori': 45, 'thi': 0.5},
    #           {'mat': m1, 'ori': 90, 'thi': 0.5},
    #           {'mat': m1, 'ori': 0, 'thi': 0.5}]
    #
    # illustrateLayup(layup1, size=(4, 3))


def illustrateCurvatures(Kx, Ky, Kxy):
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from matplotlib.ticker import LinearLocator, FormatStrFormatter
    import numpy as np
    fig = plt.figure(figsize=(6, 4))
    ax = fig.gca(projection='3d')
    # Make data.
    X = np.arange(-0.5, 0.6, 0.1)
    Y = np.arange(-0.5, 0.6, 0.1)
    X, Y = np.meshgrid(X, Y)
    Z1 = (- Kx * X ** 2 - Ky * Y ** 2 - Kxy * X * Y) / 2
    surf1 = ax.plot_surface(X, Y, Z1, cmap=cm.coolwarm,
                            linewidth=0, antialiased=False)

    # # EXAMPLE:
    # % matplotlib
    # inline
    # illustrateCurvatures(Kx=-0.5, Ky=0.5, Kxy=0.0)
    #


def plotLayerStresses(results):
    h = []
    sx, sy, sxy = [], [], []
    s1, s2, s12 = [], [], []
    for layer in results:
        h.append(layer['h']['bot'])
        h.append(layer['h']['top'])
        sx.append(layer['stress']['xyz']['bot'][0])
        sx.append(layer['stress']['xyz']['top'][0])
        sy.append(layer['stress']['xyz']['bot'][1])
        sy.append(layer['stress']['xyz']['top'][1])
        sxy.append(layer['stress']['xyz']['bot'][2])
        sxy.append(layer['stress']['xyz']['top'][2])
        s1.append(layer['stress']['123']['bot'][0])
        s1.append(layer['stress']['123']['top'][0])
        s2.append(layer['stress']['123']['bot'][1])
        s2.append(layer['stress']['123']['top'][1])
        s12.append(layer['stress']['123']['bot'][2])
        s12.append(layer['stress']['123']['top'][2])
    import matplotlib.pyplot as plt
    fig, (ax1, ax2) = plt.subplots(ncols=2, nrows=1, figsize=(10, 4))
    ax1.grid(True)
    ax1.plot(sx, h, '-', color='red', label='$\sigma_x$')
    ax1.plot(sy, h, '-', color='blue', label='$\sigma_y$')
    ax1.plot(sxy, h, '-', color='green', label=r'$\tau_{xy}$')
    ax1.set_xlabel('Stress', fontsize=12)
    ax1.set_ylabel('z', fontsize=14)
    ax1.legend(loc='best')
    ax2.grid(True)
    ax2.plot(s1, h, '--', color='red', label='$\sigma_1$')
    ax2.plot(s2, h, '--', color='blue', label='$\sigma_2$')
    ax2.plot(s12, h, '--', color='green', label=r'$\tau_{12}$')
    ax2.set_xlabel('Stress', fontsize=12)
    ax2.set_ylabel('z', fontsize=14)
    ax2.legend(loc='best')
    plt.tight_layout()
    plt.show()

    #
    # # EXAMPLE:
    # from laminatelib import laminateStiffnessMatrix, solveLaminateLoadCase, \
    #     layerResults
    #
    # m1 = {'E1': 140000, 'E2': 10000, 'v12': 0.3, 'G12': 5000, 'XT': 1200, 'XC': 800,
    #       'YT': 50, 'YC': 120, 'S12': 75, 'S23': 50, 'f12': -0.5}
    #
    # layup1 = [{'mat': m1, 'ori': 0, 'thi': 1},
    #           {'mat': m1, 'ori': 90, 'thi': 1},
    #           {'mat': m1, 'ori': 0, 'thi': 1}]
    # ABD1 = laminateStiffnessMatrix(layup1)
    # load, defs = solveLaminateLoadCase(ABD=ABD1, Nx=500, Nxy=200)
    # res = layerResults(layup1, defs)
    # % matplotlib
    # inline
    # plotLayerStresses(res)
    #


def plotLayerFailure(results):
    h = []
    ms, me, tw = [], [], []
    for layer in results:
        h.append(layer['h']['bot'])
        h.append(layer['h']['top'])
        ms.append(layer['fail']['MS']['bot'])
        ms.append(layer['fail']['MS']['top'])
        me.append(layer['fail']['ME']['bot'])
        me.append(layer['fail']['ME']['top'])
        tw.append(layer['fail']['TW']['bot'])
        tw.append(layer['fail']['TW']['top'])
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(5, 4))
    ax.grid(True)
    ax.plot(ms, h, '-', color='red', label='$f_E (MS)$')
    ax.plot(me, h, '-', color='blue', label='$f_E (ME)$')
    ax.plot(tw, h, '-', color='green', label='$f_E (TW)$')
    ax.set_xlabel('$f_E$', fontsize=12)
    ax.set_ylabel('z', fontsize=14)
    ax.legend(loc='best')
    plt.tight_layout()
    plt.show()


# EXAMPLE:
# from laminatelib import laminateStiffnessMatrix, solveLaminateLoadCase, \
#    layerResults

# m1 = {'E1': 140000, 'E2': 10000, 'v12': 0.3, 'G12': 5000, 'XT': 1200, 'XC': 800,
#       'YT': 50, 'YC': 120, 'S12': 75, 'S23': 50, 'f12': -0.5}
#
# layup1 = [{'mat': m1, 'ori': 0, 'thi': 1},
#           {'mat': m1, 'ori': 90, 'thi': 1},
#           {'mat': m1, 'ori': 0, 'thi': 1},
#           {'mat': m1, 'ori': 90, 'thi': 1},
#           {'mat': m1, 'ori': 0, 'thi': 1}]

    # ABD1 = laminateStiffnessMatrix(layup1)
    # load, defs = solveLaminateLoadCase(ABD=ABD1, Ny=500, Mx=200)
    # res = layerResults(layup1, defs)
    # % matplotlib
    # inline
    # plotLayerFailure(res)

def plotHashin(results):
    h = []
    ff, iff = [], []
    for layer in results:
        h.append(layer['h']['bot'])
        h.append(layer['h']['top'])
        ff.append(layer['fail']['HNFF']['bot'])
        ff.append(layer['fail']['HNFF']['top'])
        iff.append(layer['fail']['HNIFF']['bot'])
        iff.append(layer['fail']['HNIFF']['top'])
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(5, 4))
    ax.grid(True)
    ax.plot(ff, h, '-', color='red', label='$f_E (ff)$')
    ax.plot(iff, h, '-', color='blue', label='$f_E (iff)$')
    ax.set_xlabel('$f_E$', fontsize=12)
    ax.set_ylabel('z', fontsize=14)
    ax.legend(loc='best')
    plt.tight_layout()
    plt.show()


# EXAMPLE:
#from laminatelib import laminateStiffnessMatrix, solveLaminateLoadCase, \
#    layerResults

m1 = {'E1': 140000, 'E2': 10000, 'v12': 0.3, 'G12': 5000, 'XT': 1200, 'XC': 800,
      'YT': 50, 'YC': 120, 'S12': 75, 'S23': 50, 'f12': -0.5}

layup1 = [{'mat': m1, 'ori': 0, 'thi': 1},
          {'mat': m1, 'ori': 90, 'thi': 1},
          {'mat': m1, 'ori': 0, 'thi': 1},
          {'mat': m1, 'ori': 90, 'thi': 1},
          {'mat': m1, 'ori': 0, 'thi': 1}]

    # ABD1 = laminateStiffnessMatrix(layup1)
    # load, defs = solveLaminateLoadCase(ABD=ABD1, Ny=500, Mx=200)
    # res = layerResults(layup1, defs)
    # % matplotlib
    # inline
    # plotHashin(res)

# import numpy as np
# import matplotlib.pyplot as plt
# %matplotlib inline
# from scipy.stats import norm
#
# fig,ax = plt.subplots(figsize=(8,4))
#
# meanL  = 15
# covL   = 0.1
# stdL   = covL*meanL
# xL=np.linspace(5,25,1000)
# fL=norm.pdf(xL, meanL, stdL)
# ax.plot(xL, fL, color='blue', label='Load')
#
# ax.set_xlabel('Load (kg)')
# ax.set_ylabel('Frequency')
# ax.set_xlim(0,)
# ax.set_ylim(0,)
# ax.grid(True)
# ax.legend(loc='best')
#
# meanR  = 18
# covR   = 0.08
# stdR   = covR*meanR
# xR=np.linspace(5,30,1000)
# fR=norm.pdf(xR, meanR, stdR)
# ax.plot(xR, fR, color='red', label='Resistance')
#
# ax.set_xlabel('Load, Resistance (kg)')
# ax.set_ylabel('Frequency')
# ax.set_xlim(0,)
# ax.set_ylim(0,0.35)
# ax.grid(True)
# ax.legend(loc='best')
#
# plt.tight_layout()
# plt.savefig('imgs/safety-4.png')
# plt.show()