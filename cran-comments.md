## Third Submission

* Incorporated Comment (Please do not modifiy the .GlobalEnv. This is not allowed by the CRAN
policies. -> R/load_shapeaustria.R) by removing the function and adapting and reorganizing the example.

### Test environments
* Ubuntu 20.04.4 LTS, R 4.2.0
* Fedora Linux, R-devel, clang, gfortran

### R CMD check results (locally)
0 errors ✔ | 0 warnings ✔ | 0 notes ✔


## Second Resubmission

* Removed examples with CPU time > 2.5 times elapsed time and added automated tests instead.

### Test environments
* Ubuntu 20.04.4 LTS, R 4.2.0
* Fedora Linux, R-devel, clang, gfortran

### R CMD check results (locally)
0 errors ✔ | 0 warnings ✔ | 0 notes ✔

## Second Submission

Thank you for the valuable comments. I did my best to incorporate your comments

* __Comment1__: Please always write package names, software names ect. in single quotes in title and description.
__A__: checked and adpated accordingly

* __Comment2__: References in the Description File
__A__: checked and adpated accordingly

* __Comment3__: Write TRUE and FALSE instead of T and F.
__A__: checked and adpated accordingly

* __Comment4__: Missing value to .Rd for 
      print.SAEforest.Rd: value
      tune_parameters.Rd: value
__A__: Added a value term

* __Comment5__: dontrun{} and examples
__A__: replaced dontrun{} with donttest whenever compuationally reasonable. 
Additionally modified some examples such that automated testing is possible in < 5 sec.

* __Comment6__: Information messages to the console:
__A__: Checked and adapted.

* __Comment7__: Please do not modifiy the .GlobalEnv. This is not allowed by the CRAN
policies. -> R/load_shapeaustria.R
__A__: There is no modifaction of .GlobalEnv in the code. This function loads a shape-file for
 an example for users to illustrate and plot results on a map. Function is not needed
 per se for the package. 

### Test environments
* Ubuntu 20.04.4 LTS, R 4.2.0
* Fedora Linux, R-devel, clang, gfortran

### R CMD check results (locally)
0 errors ✔ | 0 warnings ✔ | 0 notes ✔


## First Submission

### Test environments
* Ubuntu 20.04.4 LTS, R 4.2.0
* Fedora Linux, R-devel, clang, gfortran

### R CMD check results (locally)
0 errors ✔ | 0 warnings ✔ | 0 notes ✔
