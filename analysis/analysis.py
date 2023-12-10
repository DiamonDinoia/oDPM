import pyodpm
import numpy as np
import matplotlib.pyplot as plt


histories = int(10e5)


def convertArray(geometry):
    array = pyodpm.toNdArray(geometry)
    return np.ravel(array).reshape(array.shape[0], array.shape[1], array.shape[2], order='F')


def getLateralProfile(array):
    x, y = int(array.shape[0]/2), int(array.shape[1]/2)
    return list(array[x, y, :])


def getFrontalProfile(array):
    y, z = int(array.shape[1]/2), int(array.shape[2]/2)
    return list(array[:, y, z])


def getSideProfile(array):
    x, z = int(array.shape[0]/2), int(array.shape[2]/2)
    return list(array[x, :, z])


def getDepthDistribution(array):
    return np.sum(np.sum(array, axis=0, initial=0), axis=0, initial=0)


run = pyodpm.G4Run([1, 2], pyodpm.G4DATASET_DIR)
run.initializeWater(33, 33, 33, 1)
run.initializeGeometry()
reference = run.Run("e-", pyodpm.ThreeVector_d(0.0, 0.0, 0.001),
                    pyodpm.ThreeVector_d(0., 0., 1.), 6*pyodpm.MeV, histories)
geometry = run.convertGeometry()
dpm_run = pyodpm.Run(42)
beam = pyodpm.ElectronPencilBeam(pyodpm.ThreeVector_d(
    0.0, 0.0, -1.*0.49), pyodpm.ThreeVector_d(0., 0., 1.), 6*pyodpm.MeV)
pyodpm.simulate(dpm_run, run, geometry, beam, histories)
maxError, avgError, gammaPassingRate, minGamma = pyodpm.validate(reference, geometry)

print(maxError, avgError, gammaPassingRate, minGamma)

reference_array = convertArray(reference)
geometry_array = convertArray(geometry)


plt.plot(getDepthDistribution(reference_array), label="Geant4")
plt.plot(getDepthDistribution(geometry_array), label="ODPM")
plt.legend()
plt.show()
plt.savefig("depth.png")

plt.plot(getLateralProfile(reference_array), label="Geant4")
plt.plot(getLateralProfile(geometry_array), label="ODPM")
plt.legend()
plt.show()
plt.savefig("lateral.png")

plt.plot(getFrontalProfile(reference_array), label="Geant4")
plt.plot(getFrontalProfile(geometry_array), label="ODPM")
plt.legend()
plt.show()
plt.savefig("frontal.png")

plt.plot(getSideProfile(reference_array), label="Geant4")
plt.plot(getSideProfile(geometry_array), label="ODPM")
plt.legend()
plt.show()
plt.savefig("side.png")
