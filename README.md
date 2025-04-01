# 0-1D cardiovascular solver for the publication "Personalized Pulse Wave Propagation Modeling to Improve Vasopressor Dosing Management in Patients with Severe Traumatic Brain Injury"

## Installation
The presentend solver is written in C++. The wrapper is written in Python (3.13.0).

1. Clone repository.
2. Create python environment.
3. Change python environment to already created. Install requirements from the file `requirements.txt` using command `pip install -r requirements.txt`
4. Go to the file PythonPWAExtension_library → SolverSpeed → setup.py. It is the installation file with python wrapper to cpp solver.
5. Open setup.py with editor and change extra_compile_args and extra_link_args to paths with gsl library (if needed).
6. Install solver as the python library using command `pip install .`.
7. After installation you should receive following message: `Successfully installed PythonPWAExtension-0.9.8`
   - to be sure, you can check if library is installed in the given env. To do this run (in the given env!) python3 and execute following commands:
    ```
    >>> import pkg_resources
    >>> pkg_resources.require("PythonPWAExtension")[0].version
    ```
    You should receive message: `0.9.8`.
 8. Additionaly you can change number of threads in `omp_set_num_threads` in the file `PythonPWAExtension_library` &rarr; `SolverSpeed` &rarr; `C++ source code` &rarr; `PythonPWAspeed.cpp` (if needed).

 ## Sensitive data
 To get the data patients, please write an e-mail to the authors of the publication (Kamil Wolos or Jan Poleszczuk).