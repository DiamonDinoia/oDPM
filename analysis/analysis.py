import matplotlib.pyplot as plt
import numpy as np

import pyodpm

histories = int(1_000_00)

def convertArray(geometry):
    array = pyodpm.toNdArray(geometry)
    return np.ravel(array).reshape(array.shape[0], array.shape[1], array.shape[2], order='F')


def getLateralProfile(array):
    x, y = int(array.shape[0] / 2), int(array.shape[1] / 2)
    return list(array[x, y, :])


def getFrontalProfile(array):
    y, z = int(array.shape[1] / 2), int(array.shape[2] / 2)
    return list(array[:, y, z])


def getSideProfile(array):
    x, z = int(array.shape[0] / 2), int(array.shape[2] / 2)
    return list(array[x, :, z])

def getDepthDistribution(array):
    return np.sum(np.sum(array, axis=0, initial=0), axis=0, initial=0)


run = pyodpm.G4Run([1, 2], pyodpm.G4DATASET_DIR, 'FTFP_BERT_Fast')
run.initializeWater(333, 333, 333, 1)
# run.initializeDicom("/home/mbarbone/repos/photonMonteCarlo/lung_ct/1-", 142)
center = run.center()
# print("center=", center.x, center.y, center.z)
# if center is 0 the next code will segfault
reference = run.Run("gamma", pyodpm.ThreeVector_d(0.0, 0.0, -center.z),
                    pyodpm.ThreeVector_d(0., 0., 1.), 6 * pyodpm.MeV, histories)

geometry = run.convertGeometry()

dpm_run = pyodpm.Run(42)
beam = pyodpm.PhotonPencilBeam(pyodpm.ThreeVector_d(center.x, center.y, 0),
                                 pyodpm.ThreeVector_d(0., 0., 1.), 6 * pyodpm.MeV)
pyodpm.simulate(dpm_run, run, geometry, beam, histories)
# maxError, avgError, gammaPassingRate, minGamma = pyodpm.validate(reference, geometry)

# print(maxError, avgError, gammaPassingRate, minGamma)

reference_array = convertArray(reference)
geometry_array = convertArray(geometry)

plt.plot(getDepthDistribution(reference_array), label="Geant4")
plt.plot(getDepthDistribution(geometry_array), label="ODPM")
plt.legend()
plt.savefig("e_depth.png")
plt.draw()
plt.show()
plt.clf()

plt.plot(getLateralProfile(reference_array), label="Geant4")
plt.plot(getLateralProfile(geometry_array), label="ODPM")
plt.legend()
plt.savefig("e_lateral.png")
plt.draw()
plt.show()
plt.clf()

plt.plot(getFrontalProfile(reference_array), label="Geant4")
plt.plot(getFrontalProfile(geometry_array), label="ODPM")
plt.legend()
plt.savefig("e_frontal.png")
plt.draw()
plt.show()
plt.clf()

plt.plot(getSideProfile(reference_array), label="Geant4")
plt.plot(getSideProfile(geometry_array), label="ODPM")
plt.legend()
plt.savefig("e_side.png")
plt.draw()
plt.show()
plt.clf()

plt.imshow(geometry_array[:, :, 0], label="ODPM", cmap='Blues')
plt.draw()
plt.show()
plt.clf()
plt.imshow(reference_array[:, :, 1], label="Geant4", cmap='Reds')
plt.draw()
plt.show()
plt.legend()
plt.savefig("e_scatter.png")

plt.clf()
