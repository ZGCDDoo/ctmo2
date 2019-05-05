.. _params-example:

Example Parameter file
=======================

.. code-block:: json

    {
        "ctmoVersion": "2.0.0",
        "logging": {
            "file": "ctmo.log",
            "level": "debug",
            "logToFile": true
        },
        "model": {
            "J_H": 0.0,
            "U": 20.0,
            "UPrime": 0.0,
            "beta": 5.0,
            "cluster": {
                "Nx": 1,
                "Ny": 1,
                "Nz": 1
            },
            "delta": 0.01,
            "gPhonon": 0.0,
            "hybUpFile": "hybUp50.dat",
            "modelFile": "SIAM_Square.model",
            "mu": 10.0,
            "nOrb": 1,
            "nkpts": 256,
            "tParameters": {
                "00": {
                    "t2x": 0.0,
                    "t2y": 0.0,
                    "t2z": 0.0,
                    "t3": 0.0,
                    "tIntra": 0.0,
                    "tx": -1.0,
                    "tx=-y": 0.0,
                    "tx=-z": 0.0,
                    "tx=y": 0.0,
                    "tx=z": 0.0,
                    "ty": -1.0,
                    "ty=-z": 0.0,
                    "ty=z": 0.0,
                    "tz": 0.0
                }
            },
            "w0Phonon": 0.0
        },
        "monteCarlo": {
            "measurementTime": 30.0,
            "seed": 203937213,
            "thermFromConfig": false,
            "thermalizationTime": 0.5
        },
        "selfCon": {
            "eCutSelfCon": 200,
            "weightsI": 0.2,
            "weightsR": 0.2
        },
        "slmc": {
            "batchSize": 500,
            "maxBatchesSaved": 40,
            "measurementTime": 600.0,
            "simulationId": -9999,
            "thermalizationTime": 5.0,
            "updatesMeas": 50000
        },
        "solver": {
            "S": 1.0,
            "averageOrbitals": true,
            "cleanUpdate": 10000,
            "eCutGreen": 100,
            "isOneOrbitalOptimized": true,
            "isOrbitalDiagonal": true,
            "n_tau_sampling": 5,
            "ntau": 5000,
            "updatesMeas": 100
        }
    }
