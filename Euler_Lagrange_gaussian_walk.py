from random import gauss
from math import sqrt
from numpy import identity,matrix
from numpy.linalg import cholesky,LinAlgError
import scipy.io


def LTsolve(a,b):
    d = len(a)
    x = matrix([[0.0]]*d)
    for i in range(d):
        x[i,0] = (b[i,0]-((a[i,:i]*x[:i])[0,0]))/(a[i,i])
    return x

def eulerlagrangefield(point,tangent):
    mu,sigma = point
    dmu,dsigma = tangent
    logarithmicderivative = dsigma*(sigma.getI())
    return (logarithmicderivative*dmu,(logarithmicderivative*dsigma)-(dmu*(dmu.getT())))

boundaryhits = 0

def cautiouslysimulatethendecompose(initialvalue,initialtangent,accelfield,stepsize,numsteps):
    # This uses Euler's method with ordered pairs of matrices.
    currentvalue = initialvalue
    currenttangent = initialtangent
    nextvalue = [currentvalue[i]+(stepsize*currenttangent[i]) for i in range(2)]
    try:
        oldesky = cholesky(nextvalue[1])
    except LinAlgError:
        print('reached boundary')
        boundaryhits += 1
        return (nextvalue[0],matrix(identity(len(initialtangent[0]))))
    for t in range(numsteps-2):
        # the main loop of Euler's method
        doublederiv = accelfield(currentvalue,currenttangent)
        currentvalue = nextvalue
        currenttangent = [currenttangent[i]+(stepsize*doublederiv[i]) for i in range(2)]
        nextvalue = [currentvalue[i]+(stepsize*currenttangent[i]) for i in range(2)]
        # the exception-handling is the easiest way to quickly make
        try:
            temp = cholesky(nextvalue[1])
            oldesky = temp
        except LinAlgError:
            print('reached boundary')
            boundaryhits += 1
            return (nextvalue[0],oldesky)
        # sure the "covariance" matrices are staying positive-definite
    doublederiv = accelfield(currentvalue,currenttangent)
    currentvalue = nextvalue
    currenttangent = [currenttangent[i]+(stepsize*doublederiv[i]) for i in range(2)]
    nextvalue = [currentvalue[i]+(stepsize*currenttangent[i]) for i in range(2)]
    try:
        temp = cholesky(nextvalue[1])
        return [nextvalue[0],temp]
    except LinAlgError:
        print('reached boundary')
        boundaryhits += 1
        return (nextvalue[0],oldesky)
    # That try-except block immediately above this line had better cause this function to return.
    # If the next line somehow executes, then I have no clue what went wrong.
    raise Exception

def simulatethendecompose(initialvalue,initialtangent,accelfield,stepsize,numsteps):
    # This uses Euler's method with ordered pairs of matrices.
    currentvalue = initialvalue
    currenttangent = initialtangent
    doublederiv = accelfield(currentvalue,currenttangent)
    currentvalue = [currentvalue[i]+(stepsize*currenttangent[i]) for i in range(2)]
    for t in range(numsteps-1):
        # the main loop of Euler's method
        currenttangent = [currenttangent[i]+(stepsize*doublederiv[i]) for i in range(2)]
        doublederiv = accelfield(currentvalue,currenttangent)
        currentvalue = [currentvalue[i]+(stepsize*currenttangent[i]) for i in range(2)]
    return currentvalue

def gaussianwalk(dimension,datasamplesize,outputsize,oderunsperpoint,stepsperoderun):
    if datasamplesize <= 0 or oderunsperpoint < 1 or stepsperoderun < 1:
        print('invalid datasamplesize or oderunsperpoint or stepsperoderun')
        raise(ValueError)
    if stepsperoderun < 2:
        print('I optimized this for  2 <= stepsperoderun  at the cost of making it not handle smaller values of stepsperoderun correctly.')
        raise(NotImplementedError)
    h = 1/stepsperoderun
    s = sqrt(1.0/(datasamplesize*oderunsperpoint)) # The oderuns are independent, and I think 1/datasamplesize is the variance that the points are supposed to have.
    stretch = (sqrt(2))*s # I'll be using this because the tensor at the base-point is not the identity matrix,
    zeroarray = [[0.0 for i in range(dimension)] for j in range(dimension)]
    basemu = matrix([[0.0]]*dimension)
    baseproduct = matrix(identity(dimension))
    for outputnumber in range(outputsize):
        currentmu = basemu
        bigproduct = baseproduct
        for odenumber in range(oderunsperpoint):
            mudirection = matrix([[gauss(0,s)] for _ in range(dimension)])
            sigmadirection = zeroarray
            for j in range(dimension):
                sigmadirection [j][j] = gauss(0,stretch) # and I don't want to recompute (sqrt(2))*s each time.
                for i in range(j):
                    sigmadirection[i][j] = gauss(0,s)
                    sigmadirection[j][i] = sigmadirection[i][j]
            deltamu,L = cautiouslysimulatethendecompose([basemu,baseproduct],[mudirection,matrix(sigmadirection)],eulerlagrangefield,h,stepsperoderun)
            currentmu = L*(currentmu+LTsolve(L,deltamu))
            bigproduct = L*bigproduct
        yield [currentmu,bigproduct*(bigproduct.getT())]





def randomCov(dimension, datasamplesize, outputsize, oderunsperpoint=128, stepsperoderun=32):
    Sall = [[[0.0]*3]*3]*outputsize
    index = 0
    for point in gaussianwalk(dimension,datasamplesize,outputsize,oderunsperpoint, stepsperoderun): 
        M,S = point      
        print(M)        
        print(S)
        Sall[index] = S
        index += 1
    print(str(boundaryhits)+' boundary hits')
    filename = str(dimension) + ',' + str(datasamplesize) + "," + str(outputsize) + "," + str(oderunsperpoint) + "," + str(stepsperoderun) + ".mat"
    return (Sall, filename)


# Run code and save output to file
outputs, filename = randomCov(3,4,1000)
scipy.io.savemat('/Users/drwho/Desktop/Ricky/'+filename, mdict={'cov': outputs})

