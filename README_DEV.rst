==========================================================================
 CTMO2.0 : Continuous time Multiorbital
==========================================================================

:Authors: Charles-David Hébert 
:Date: $Date: 2018-02-06 $
:Revision: $Revision: 2.0.0 $
:Description: Description

Naming conventions
-------------------
The naming conventions and style is mainly derived from the Visual studio style for C++.
See this site:
    https://docs.microsoft.com/en-us/dotnet/standard/design-guidelines/naming-guidelines
 

General guidelines before pull requests
----------------------------------------

Build
^^^^^^^^^^^^^^^^^^^^^^
* Must build with clang3.8-clang6.0 for serial mode, no warnings, except for external libraries
* Must build with g++5.4-g++8.0 for MPI mode, no warings.
* Must pass all unit tests
* Must pass the integration tests


Naming conventions
^^^^^^^^^^^^^^^^^^^
* Must respect the naming conventions


Formatting
^^^^^^^^^^^^^^^^
* Formatted according to Visual Studio standard, for exemple in visual studio code, put the folling preferences:
* "C_Cpp.clang_format_style": "Visual Studio"
* "editor.formatOnSave": true
* Or run the script format.sh: bash format.sh
* Check the warnings that the script cppcheck.sh outputs.



    