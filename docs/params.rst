.. _params:


The Parameter file
===================

The parameter files used to specify the details of the Monte Carlo simulations are in the json format. 
This file format is convinient, readable and also ubiquitous. While reading the following descriptions, please open up the file 
‘test/Simulations/holstein_DMFT_T0.3/params1.json’ to follow along. The name of the parameters follow the simple camelCase convention.

You can also reference the parameter file here: :ref:`params-example`


monteCarlo
-----------

The parameters in this sub-json specify the monte carlo specific details:

.. describe:: measurementTime (float)

        The time in minutes each processor measures.


.. describe:: thermalizationTime (float)

        The time in minutes for which each processor will thermalize. 
        It is difficult to give a good optimal value. I would say, ~10% of the measurement time.

.. describe:: seed (int)

        The seed for the random number generator.

.. describe:: thermalizeFromConfig (bool)

        Should the program thermalize or simply start measuring format the previously saved configuraion. 
        If true thermalizes from the previously saved configuration.


model
-----



 