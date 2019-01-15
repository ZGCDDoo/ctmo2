.. _params-example:

Example Parameter file
=======================

.. code-block:: json

    {
        "monteCarlo": {
            "measurementTime": 5.0,
            "thermalizationTime": 0.5,
            "thermFromConfig": false,
            "seed": 1096954689
        },
        "model": {
            "beta": 33.3333333,
            "hybUpFile": "hybUp1.dat",
            "modelFile": "SIAM_Square.model",
            "J_H": 0.0,
            "nkpts": 60,
            "nOrb": 1,
            "U": 0.0,
            "UPrime": 0.0,
            "delta": 0.01,
            "gPhonon": 0.08,
            "w0Phonon": 0.08,
            "mu": 0.0,
            "cluster": {
                "Nx": 1,
                "Ny": 1,
                "Nz": 1
            },
            "tParameters": {
                "00": {
                    "t2x": 0.0,
                    "t2y": 0.0,
                    "t2z": 0.0,
                    "tIntra": 0.0,
                    "tx": -0.3926,
                    "tx=-y": 0.0,
                    "tx=-z": 0.0,
                    "tx=y": 0.0,
                    "tx=z": 0.0,
                    "ty": -0.3926,
                    "ty=-z": 0.0,
                    "ty=z": 0.0,
                    "tz": -0.3926,
                    "t3": 0.0
                }
            }
        },
        "solver": {
            "cleanUpdate": 100000,
            "updatesMeas": 50,
            "averageOrbitals": true,
            "isOrbitalDiagonal": true,
            "eCutGreen": 100,
            "ntau": 5000,
            "n_tau_sampling": 5
        },
        "selfCon": {
            "eCutSelfCon": 200,
            "weightsR": 0.2,
            "weightsI": 0.2
        },
        "logging": {
            "level": "debug",
            "logToFile": false,
            "file": "ctmo.log"
        }
    }