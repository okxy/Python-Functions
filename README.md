#Python Functions

A systematically curated list of Python functions, modules, frameworks. Enviouly inspired by [Matlab Documentation Center](http://www.mathworks.com/help/matlab/index.html)


- [Getting Started with Python](#Getting-Started-with-Python)
    - [Environment Management](#environment-management)
    - [Package Management](#package-management)
   
    - [Miscellaneous](#miscellaneous)
    - [Algorithms and Design Patterns](#algorithms-and-design-patterns)
    - [Editor Plugins](#editor-plugins)
- [Resources](#resources)
    - [Websites](#websites)
    - [Weekly](#weekly)
    - [Twitter](#twitter)
- [Other Awesome Lists](#other-awesome-lists)
- [Contributing](#contributing)

---

## Getting Started with Python

*Basic functions about Python.*

* [help](help) - help() function
* [dir](builtin) - function
* 

## Language Fundamentals


## Mathematics

###Elementary Math
###Linear Algebra
###Statistics and Random Numbers
###Interpolation
###Optimization
###Numerical Integration and Differential Equation
###Fourier Analysis and Filtering
###Sparse Matrices
###COmputational Geometry
## Graphics
## Programming Scripts and Functions
## Data and File Management
## GUI Building
## Advanced Software Development
## Desktop Environment

## Environment Management

*Libraries for Python version and environment management.*

* [pyenv](https://github.com/yyuu/pyenv) - Simple Python version management.
* [virtualenv](https://pypi.python.org/pypi/virtualenv) - A tool to create isolated Python environments.
* [virtualenvwrapper](https://pypi.python.org/pypi/virtualenvwrapper) - A set of extensions to virtualenv
* [virtualenv-api](https://github.com/sjkingo/virtualenv-api) - An API for virtualenv and pip.
* [pew](https://pypi.python.org/pypi/pew/) - A set of tools to manage multiple virtual environments.
* [Vex](https://github.com/sashahart/vex) - Run a command in the named virtualenv.
* [PyRun](https://www.egenix.com/products/python/PyRun/) - A one-file, no-installation-needed version of Python.

## Package Management

Language Fundamentals

Entering Commands
ans	Most recent answer
clc	Clear Command Window
diary	Save Command Window text to file
format	Set display format for output
home	Send cursor home
iskeyword	Determine whether input is MATLAB keyword
more	Control paged output for Command Window
Matrices and Arrays
Array Creation and Concatenation

accumarray	Construct array with accumulation
blkdiag	Construct block diagonal matrix from input arguments
diag	Get diagonal elements or create diagonal matrix
eye	Identity matrix
false	Logical 0 (false)
freqspace	Frequency spacing for frequency response
linspace	Generate linearly spaced vectors
logspace	Generate logarithmically spaced vectors
meshgrid	Rectangular grid in 2-D and 3-D space
ndgrid	Rectangular grid in N-D space
ones	Create array of all ones
rand	Uniformly distributed pseudorandom numbers
true	Logical 1 (true)
zeros	Create array of all zeros
cat	Concatenate arrays along specified dimension
horzcat	Concatenate arrays horizontally
vertcat	Concatenate arrays vertically
Indexing

colon	Create vectors, array subscripting, and for-loop iterators
end	Terminate block of code, or indicate last array index
ind2sub	Subscripts from linear index
sub2ind	Convert subscripts to linear indices
Array Dimensions

length	Length of vector or largest array dimension
ndims	Number of array dimensions
numel	Number of array elements
size	Array dimensions
height	Number of table rows
width	Number of table variables
iscolumn	Determine whether input is column vector
isempty	Determine whether array is empty
ismatrix	Determine whether input is matrix
isrow	Determine whether input is row vector
isscalar	Determine whether input is scalar
isvector	Determine whether input is vector
Sorting and Reshaping Arrays

blkdiag	Construct block diagonal matrix from input arguments
circshift	Shift array circularly
ctranspose	Complex conjugate transpose
diag	Get diagonal elements or create diagonal matrix
flip	Flip order of elements
fliplr	Flip array left to right
flipud	Flip array up to down
ipermute	Inverse permute dimensions of N-D array
permute	Rearrange dimensions of N-D array
repmat	Replicate and tile array
reshape	Reshape array
rot90	Rotate array 90 degrees
shiftdim	Shift dimensions
issorted	Determine whether set elements are in sorted order
sort	Sort array elements
sortrows	Sort array rows
squeeze	Remove singleton dimensions
transpose	Transpose
vectorize	Vectorize expression
Operators and Elementary Operations
Arithmetic

plus	Addition
uplus	Unary plus
minus	Subtraction
uminus	Unary minus
times	Element-wise multiplication
rdivide	Right array division
ldivide	Left array division
power	Element-wise power
mtimes	Matrix Multiplication
mrdivide	Solve systems of linear equations xA = B for x
mldivide	Solve systems of linear equations Ax = B for x
mpower	Matrix power
cumprod	Cumulative product
cumsum	Cumulative sum
diff	Differences and Approximate Derivatives
prod	Product of array elements
sum	Sum of array elements
ceil	Round toward positive infinity
fix	Round toward zero
floor	Round toward negative infinity
idivide	Integer division with rounding option
mod	Modulus after division
rem	Remainder after division
round	Round to nearest integer
Relational Operations

Relational Operators	Relational operations
eq	Determine equality
ge	Determine greater than or equal to
gt	Determine greater than
le	Determine less than or equal to
lt	Determine less than
ne	Determine inequality
isequal	Determine array equality
isequaln	Determine array equality, treating NaN values as equal
Logical Operations

Logical Operators: Short-circuit	Logical operations with short-circuiting
and	Find logical AND
not	Find logical NOT
or	Find logical OR
xor	Logical exclusive-OR
all	Determine if all array elements are nonzero or true
any	Determine if any array elements are nonzero
false	Logical 0 (false)
find	Find indices and values of nonzero elements
islogical	Determine if input is logical array
logical	Convert numeric values to logicals
true	Logical 1 (true)
Set Operations

intersect	Set intersection of two arrays
ismember	Array elements that are members of set array
issorted	Determine whether set elements are in sorted order
setdiff	Set difference of two arrays
setxor	Set exclusive OR of two arrays
union	Set union of two arrays
unique	Unique values in array
join	Merge two tables by matching up rows using key variables
innerjoin	Inner join between two tables
outerjoin	Outer join between two tables
Bit-Wise Operations

bitand	Bit-wise AND
bitcmp	Bit-wise complement
bitget	Get bit at specified position
bitor	Bit-wise OR
bitset	Set bit at specific location
bitshift	Shift bits specified number of places
bitxor	Bit-wise XOR
swapbytes	Swap byte ordering
Special Characters
Special Characters	Special characters
colon	Create vectors, array subscripting, and for-loop iterators
Data Types
Numeric Types

double	Convert to double precision
single	Convert to single precision
int8	Convert to 8-bit signed integer
int16	Convert to 16-bit signed integer
int32	Convert to 32-bit signed integer
int64	Convert to 64-bit signed integer
uint8	Convert to 8-bit unsigned integer
uint16	Convert to 16-bit unsigned integer
uint32	Convert to 32-bit unsigned integer
uint64	Convert to 64-bit unsigned integer
cast	Cast variable to different data type
typecast	Convert data types without changing underlying data
isinteger	Determine if input is integer array
isfloat	Determine if input is floating-point array
isnumeric	Determine if input is numeric array
isreal	Determine if array is real
isfinite	Array elements that are finite
isinf	Array elements that are infinite
isnan	Array elements that are NaN
eps	Floating-point relative accuracy
flintmax	Largest consecutive integer in floating-point format
Inf	Infinity
intmax	Largest value of specified integer type
intmin	Smallest value of specified integer type
NaN	Not-a-Number
realmax	Largest positive floating-point number
realmin	Smallest positive normalized floating-point number
Characters and Strings

Create and Concatenate Strings
blanks	Create string of blank characters
cellstr	Create cell array of strings from character array
char	Convert to character array (string)
iscellstr	Determine whether input is cell array of strings
ischar	Determine whether item is character array
sprintf	Format data into string
strcat	Concatenate strings horizontally
strjoin	Join strings in cell array into single string
Parse Strings
ischar	Determine whether item is character array
isletter	Array elements that are alphabetic letters
isspace	Array elements that are space characters
isstrprop	Determine whether string is of specified category
sscanf	Read formatted data from string
strfind	Find one string within another
strrep	Find and replace substring
strsplit	Split string at specified delimiter
strtok	Selected parts of string
validatestring	Check validity of text string
symvar	Determine symbolic variables in expression
regexp	Match regular expression (case sensitive)
regexpi	Match regular expression (case insensitive)
regexprep	Replace string using regular expression
regexptranslate	Translate string into regular expression
Compare Strings
strcmp	Compare strings with case sensitivity
strcmpi	Compare strings (case insensitive)
strncmp	Compare first n characters of strings (case sensitive)
strncmpi	Compare first n characters of strings (case insensitive)
Change String Case, Blanks, and Justification
blanks	Create string of blank characters
deblank	Strip trailing blanks from end of string
strtrim	Remove leading and trailing white space from string
lower	Convert string to lowercase
upper	Convert string to uppercase
strjust	Justify character array
Categorical Arrays

categorical	Create categorical array
iscategorical	Determine whether input is categorical array
categories	Categories of categorical array
iscategory	Test for categorical array categories
isordinal	Determine whether input is ordinal categorical array
isprotected	Determine whether categories of categorical array are protected
addcats	Add categories to categorical array
mergecats	Merge categories in categorical array
removecats	Remove categories from categorical array
renamecats	Rename categories in categorical array
reordercats	Reorder categories in categorical array
summary	Print summary of table or categorical array
countcats	Count occurrences of categorical array elements by category
isundefined	Find undefined elements in categorical array
Tables

table	Create table from workspace variables
array2table	Convert homogeneous array to table
cell2table	Convert cell array to table
struct2table	Convert structure array to table
table2array	Convert table to homogenous array
table2cell	Convert table to cell array
table2struct	Convert table to structure array
readtable	Create table from file
writetable	Write table to file
istable	Determine whether input is table
height	Number of table rows
width	Number of table variables
summary	Print summary of table or categorical array
intersect	Set intersection of two arrays
ismember	Array elements that are members of set array
setdiff	Set difference of two arrays
setxor	Set exclusive OR of two arrays
unique	Unique values in array
union	Set union of two arrays
join	Merge two tables by matching up rows using key variables
innerjoin	Inner join between two tables
outerjoin	Outer join between two tables
sortrows	Sort array rows
stack	Stack data from multiple variables into single variable
unstack	Unstack data from single variable into multiple variables
ismissing	Find table elements with missing values
standardizeMissing	Insert missing value indicators into table
varfun	Apply function to table variables
rowfun	Apply function to table rows
Structures

struct	Create structure array
fieldnames	Field names of structure, or public fields of object
getfield	Field of structure array
isfield	Determine whether input is structure array field
isstruct	Determine whether input is structure array
orderfields	Order fields of structure array
rmfield	Remove fields from structure
setfield	Assign values to structure array field
arrayfun	Apply function to each element of array
structfun	Apply function to each field of scalar structure
table2struct	Convert table to structure array
struct2table	Convert structure array to table
cell2struct	Convert cell array to structure array
struct2cell	Convert structure to cell array
Cell Arrays

cell	Create cell array
cell2mat	Convert cell array to numeric array
cell2struct	Convert cell array to structure array
cell2table	Convert cell array to table
celldisp	Display cell array contents
cellfun	Apply function to each cell in cell array
cellplot	Graphically display structure of cell array
cellstr	Create cell array of strings from character array
iscell	Determine whether input is cell array
iscellstr	Determine whether input is cell array of strings
mat2cell	Convert array to cell array with potentially different sized cells
num2cell	Convert array to cell array with consistently sized cells
strjoin	Join strings in cell array into single string
strsplit	Split string at specified delimiter
struct2cell	Convert structure to cell array
table2cell	Convert table to cell array
Function Handles

function_handle (@)	Handle used in calling functions indirectly
feval	Evaluate function
func2str	Construct function name string from function handle
str2func	Construct function handle from function name string
localfunctions	Function handles to all local functions in MATLAB file
functions	Information about function handle
Map Containers

containers.Map	Map values to unique keys
isKey	Determine if containers.Map object contains key
keys	Identify keys of containers.Map object
remove	Remove key-value pairs from containers.Map object
values	Identify values in containers.Map object
Time Series

Time Series Basics
append	Concatenate time series objects in time dimension
get	Query timeseries object property values
getdatasamplesize	Size of data sample in timeseries object
getqualitydesc	Data quality descriptions
getsamples	Subset of time series samples using subscripted index array
plot	Plot time series
set	Set properties of timeseries object
tsdata.event	Construct event object for timeseries object
timeseries	Create timeseries object
Data Manipulation
addsample	Add data sample to timeseries object
ctranspose	Transpose timeseries object
delsample	Remove sample from timeseries object
detrend	Subtract mean or best-fit line and all NaNs from timeseries object
filter	Shape frequency content of time-series
getabstime	Extract date-string time vector into cell array
getinterpmethod	Interpolation method for timeseries object
getsampleusingtime	Extract data samples into new timeseries object
idealfilter	Apply ideal (noncausal) filter to timeseries object
resample	Select or interpolate timeseries data using new time vector
setabstime	Set times of timeseries object as date strings
setinterpmethod	Set default interpolation method for timeseries object
synchronize	Synchronize and resample two timeseries objects using common time vector
transpose	Transpose timeseries object
Event Data
addevent	Add event to timeseries object
delevent	Remove tsdata.event objects from timeseries object
gettsafteratevent	New timeseries object with samples occurring at or after event
gettsafterevent	New timeseries object with samples occurring after event
gettsatevent	New timeseries object with samples occurring at event
gettsbeforeatevent	New timeseries object with samples occurring before or at event
gettsbeforeevent	New timeseries object with samples occurring before event
gettsbetweenevents	New timeseries object with samples occurring between events
Descriptive Statistics
iqr	Interquartile range of timeseries data
max	Maximum value of timeseries data
mean	Mean value of timeseries data
median	Median value of timeseries data
min	Minimum value of timeseries data
std	Standard deviation of timeseries data
sum	Sum of timeseries data
var	Variance of timeseries data
Time Series Collections
get (tscollection)	Query tscollection object property values
isempty (tscollection)	Determine whether tscollection object is empty
length (tscollection)	Length of time vector
plot	Plot time series
set (tscollection)	Set properties of tscollection object
size (tscollection)	Size of tscollection object
tscollection	Create tscollection object
addsampletocollection	Add sample to tscollection object
addts	Add timeseries object to tscollection object
delsamplefromcollection	Remove sample from tscollection object
getabstime (tscollection)	Extract date-string time vector into cell array
getsampleusingtime (tscollection)	Extract data samples into new tscollection object
gettimeseriesnames	Cell array of names of timeseries objects in tscollection object
horzcat (tscollection)	Horizontal concatenation for tscollection objects
removets	Remove timeseries objects from tscollection object
resample (tscollection)	Select or interpolate data in tscollection using new time vector
setabstime (tscollection)	Set times of tscollection object as date strings
settimeseriesnames	Change name of timeseries object in tscollection
vertcat (tscollection)	Vertical concatenation for tscollection objects
Data Type Identification

is*	Detect state
isa	Determine if input is object of specified class
iscategorical	Determine whether input is categorical array
iscell	Determine whether input is cell array
iscellstr	Determine whether input is cell array of strings
ischar	Determine whether item is character array
isfield	Determine whether input is structure array field
isfloat	Determine if input is floating-point array
ishghandle	True for Handle Graphics object handles
isinteger	Determine if input is integer array
isjava	Determine if input is Java object
islogical	Determine if input is logical array
isnumeric	Determine if input is numeric array
isobject	Determine if input is MATLAB object
isreal	Determine if array is real
isscalar	Determine whether input is scalar
isstr	Determine whether input is character array
isstruct	Determine whether input is structure array
istable	Determine whether input is table
isvector	Determine whether input is vector
class	Determine class of object
validateattributes	Check validity of array
whos	List variables in workspace, with sizes and types
Data Type Conversion

char	Convert to character array (string)
int2str	Convert integer to string
mat2str	Convert matrix to string
num2str	Convert number to string
str2double	Convert string to double-precision value
str2num	Convert string to number
native2unicode	Convert numeric bytes to Unicode character representation
unicode2native	Convert Unicode character representation to numeric bytes
base2dec	Convert base N number string to decimal number
bin2dec	Convert binary number string to decimal number
dec2base	Convert decimal to base N number in string
dec2bin	Convert decimal to binary number in string
dec2hex	Convert decimal to hexadecimal number in string
hex2dec	Convert hexadecimal number string to decimal number
hex2num	Convert hexadecimal number string to double-precision number
num2hex	Convert singles and doubles to IEEE hexadecimal strings
table2array	Convert table to homogenous array
table2cell	Convert table to cell array
table2struct	Convert table to structure array
array2table	Convert homogeneous array to table
cell2table	Convert cell array to table
struct2table	Convert structure array to table
cell2mat	Convert cell array to numeric array
cell2struct	Convert cell array to structure array
cellstr	Create cell array of strings from character array
mat2cell	Convert array to cell array with potentially different sized cells
num2cell	Convert array to cell array with consistently sized cells
struct2cell	Convert structure to cell array
Dates and Time
datenum	Convert date and time to serial date number
datevec	Convert date and time to vector of components
datestr	Convert date and time to string format
now	Current date and time as serial date number
clock	Current date and time as date vector
date	Current date string
calendar	Calendar for specified month
eomday	Last day of month
weekday	Day of week
addtodate	Modify date number by field
etime	Time elapsed between date vectors
Mathematics

Elementary Math
Arithmetic

plus	Addition
uplus	Unary plus
minus	Subtraction
uminus	Unary minus
times	Element-wise multiplication
rdivide	Right array division
ldivide	Left array division
power	Element-wise power
mtimes	Matrix Multiplication
mrdivide	Solve systems of linear equations xA = B for x
mldivide	Solve systems of linear equations Ax = B for x
mpower	Matrix power
cumprod	Cumulative product
cumsum	Cumulative sum
diff	Differences and Approximate Derivatives
prod	Product of array elements
sum	Sum of array elements
ceil	Round toward positive infinity
fix	Round toward zero
floor	Round toward negative infinity
idivide	Integer division with rounding option
mod	Modulus after division
rem	Remainder after division
round	Round to nearest integer
Trigonometry

sin	Sine of argument in radians
sind	Sine of argument in degrees
asin	Inverse sine in radians
asind	Inverse sine in degrees
sinh	Hyperbolic sine of argument in radians
asinh	Inverse hyperbolic sine
cos	Cosine of argument in radians
cosd	Cosine of argument in degrees
acos	Inverse cosine in radians
acosd	Inverse cosine in degrees
cosh	Hyperbolic cosine
acosh	Inverse hyperbolic cosine
tan	Tangent of argument in radians
tand	Tangent of argument in degrees
atan	Inverse tangent in radians
atand	Inverse tangent in degrees
atan2	Four-quadrant inverse tangent
atan2d	Four-quadrant inverse tangent in degrees
tanh	Hyperbolic tangent
atanh	Inverse hyperbolic tangent
csc	Cosecant of input angle in radians
cscd	Cosecant of argument in degrees
acsc	Inverse cosecant in radians
acscd	Inverse cosecant in degrees
csch	Hyperbolic cosecant
acsch	Inverse hyperbolic cosecant
sec	Secant of angle in radians
secd	Secant of argument in degrees
asec	Inverse secant in radians
asecd	Inverse secant in degrees
sech	Hyperbolic secant
asech	Inverse hyperbolic secant
cot	Cotangent of angle in radians
cotd	Cotangent of argument in degrees
acot	Inverse cotangent in radians
acotd	Inverse cotangent in degrees
coth	Hyperbolic cotangent
acoth	Inverse hyperbolic cotangent
hypot	Square root of sum of squares
Exponents and Logarithms

exp	Exponential
expm1	Compute exp(x)-1 accurately for small values of x
log	Natural logarithm
log10	Common (base 10) logarithm
log1p	Compute log(1+x) accurately for small values of x
log2	Base 2 logarithm and dissect floating-point numbers into exponent and mantissa
nextpow2	Exponent of next higher power of 2
nthroot	Real nth root of real numbers
pow2	Base 2 power and scale floating-point numbers
reallog	Natural logarithm for nonnegative real arrays
realpow	Array power for real-only output
realsqrt	Square root for nonnegative real arrays
sqrt	Square root
Complex Numbers

abs	Absolute value and complex magnitude
angle	Phase angle
complex	Create complex array
conj	Complex conjugate
cplxpair	Sort complex numbers into complex conjugate pairs
i	Imaginary unit
imag	Imaginary part of complex number
isreal	Determine if array is real
j	Imaginary unit
real	Real part of complex number
sign	Signum function
unwrap	Correct phase angles to produce smoother phase plots
Discrete Math

factor	Prime factors
factorial	Factorial of input
gcd	Greatest common divisor
isprime	Determine which array elements are prime
lcm	Least common multiple
nchoosek	Binomial coefficient or all combinations
perms	All possible permutations
primes	Prime numbers less than or equal to input value
rat	Rational fraction approximation
rats	Rational output
Polynomials

poly	Polynomial with specified roots
polyder	Polynomial derivative
polyeig	Polynomial eigenvalue problem
polyfit	Polynomial curve fitting
polyint	Integrate polynomial analytically
polyval	Polynomial evaluation
polyvalm	Matrix polynomial evaluation
residue	Convert between partial fraction expansion and polynomial coefficients
roots	Polynomial roots
Special Functions

airy	Airy Functions
besselh	Bessel function of third kind (Hankel function)
besseli	Modified Bessel function of first kind
besselj	Bessel function of first kind
besselk	Modified Bessel function of second kind
bessely	Bessel function of second kind
beta	Beta function
betainc	Incomplete beta function
betaincinv	Beta inverse cumulative distribution function
betaln	Logarithm of beta function
ellipj	Jacobi elliptic functions
ellipke	Complete elliptic integrals of first and second kind
erf	Error function
erfc	Complementary error function
erfcinv	Inverse complementary error function
erfcx	Scaled complementary error function
erfinv	Inverse error function
expint	Exponential integral
gamma	Gamma function
gammainc	Incomplete gamma function
gammaincinv	Inverse incomplete gamma function
gammaln	Logarithm of gamma function
legendre	Associated Legendre functions
psi	Psi (polygamma) function
Cartesian Coordinate System Conversion

cart2pol	Transform Cartesian coordinates to polar or cylindrical
cart2sph	Transform Cartesian coordinates to spherical
pol2cart	Transform polar or cylindrical coordinates to Cartesian
sph2cart	Transform spherical coordinates to Cartesian
Constants and Test Matrices

eps	Floating-point relative accuracy
flintmax	Largest consecutive integer in floating-point format
i	Imaginary unit
j	Imaginary unit
Inf	Infinity
pi	Ratio of circle's circumference to its diameter
NaN	Not-a-Number
isfinite	Array elements that are finite
isinf	Array elements that are infinite
isnan	Array elements that are NaN
compan	Companion matrix
gallery	Test matrices
hadamard	Hadamard matrix
hankel	Hankel matrix
hilb	Hilbert matrix
invhilb	Inverse of Hilbert matrix
magic	Magic square
pascal	Pascal matrix
rosser	Classic symmetric eigenvalue test problem
toeplitz	Toeplitz matrix
vander	Vandermonde matrix
wilkinson	Wilkinson's eigenvalue test matrix
Linear Algebra
Matrix Operations

cross	Cross product
dot	Dot product
kron	Kronecker tensor product
surfnorm	Compute and display 3-D surface normals
tril	Lower triangular part of matrix
triu	Upper triangular part of matrix
transpose	Transpose
Linear Equations

cond	Condition number with respect to inversion
condest	1-norm condition number estimate
funm	Evaluate general matrix function
inv	Matrix inverse
linsolve	Solve linear system of equations
lscov	Least-squares solution in presence of known covariance
lsqnonneg	Solve nonnegative least-squares constraints problem
pinv	Moore-Penrose pseudoinverse of matrix
rcond	Reciprocal condition number
sylvester	Solve Sylvester equation AX + XB = C for X
mldivide	Solve systems of linear equations Ax = B for x
mrdivide	Solve systems of linear equations xA = B for x
Matrix Decomposition

chol	Cholesky factorization
ichol	Incomplete Cholesky factorization
cholupdate	Rank 1 update to Cholesky factorization
ilu	Sparse incomplete LU factorization
lu	LU matrix factorization
qr	Orthogonal-triangular decomposition
qrdelete	Remove column or row from QR factorization
qrinsert	Insert column or row into QR factorization
qrupdate	Rank 1 update to QR factorization
planerot	Givens plane rotation
ldl	Block LDL' factorization for Hermitian indefinite matrices
cdf2rdf	Convert complex diagonal form to real block diagonal form
rsf2csf	Convert real Schur form to complex Schur form
gsvd	Generalized singular value decomposition
svd	Singular value decomposition
Eigenvalues and Singular Values

balance	Diagonal scaling to improve eigenvalue accuracy
cdf2rdf	Convert complex diagonal form to real block diagonal form
condeig	Condition number with respect to eigenvalues
eig	Eigenvalues and eigenvectors
eigs	Largest eigenvalues and eigenvectors of matrix
gsvd	Generalized singular value decomposition
hess	Hessenberg form of matrix
ordeig	Eigenvalues of quasitriangular matrices
ordqz	Reorder eigenvalues in QZ factorization
ordschur	Reorder eigenvalues in Schur factorization
poly	Polynomial with specified roots
polyeig	Polynomial eigenvalue problem
qz	QZ factorization for generalized eigenvalues
rsf2csf	Convert real Schur form to complex Schur form
schur	Schur decomposition
sqrtm	Matrix square root
ss2tf	Convert state-space filter parameters to transfer function form
svd	Singular value decomposition
svds	Find singular values and vectors
Matrix Analysis

bandwidth	Lower and upper matrix bandwidth
cond	Condition number with respect to inversion
condeig	Condition number with respect to eigenvalues
det	Matrix determinant
isbanded	Determine if matrix is within specific bandwidth
isdiag	Determine if matrix is diagonal
ishermitian	Determine if matrix is Hermitian or skew-Hermitian
issymmetric	Determine if matrix is symmetric or skew-symmetric
istril	Determine if matrix is lower triangular
istriu	Determine if matrix is upper triangular
norm	Vector and matrix norms
normest	2-norm estimate
null	Null space
orth	Orthonormal basis for range of matrix
rank	Rank of matrix
rcond	Reciprocal condition number
rref	Reduced row echelon form
subspace	Angle between two subspaces
trace	Sum of diagonal elements
Matrix Functions

expm	Matrix exponential
logm	Matrix logarithm
sqrtm	Matrix square root
bsxfun	Apply element-by-element binary operation to two arrays with singleton expansion enabled
arrayfun	Apply function to each element of array
accumarray	Construct array with accumulation
mpower	Matrix power
Statistics and Random Numbers
Descriptive Statistics

corrcoef	Correlation coefficients
cov	Covariance matrix
max	Largest elements in array
mean	Average or mean value of array
median	Median value of array
min	Smallest elements in array
mode	Most frequent values in array
std	Standard deviation
var	Variance
Random Number Generation

rand	Uniformly distributed pseudorandom numbers
randn	Normally distributed pseudorandom numbers
randi	Uniformly distributed pseudorandom integers
randperm	Random permutation
rng	Control random number generation
RandStream	Random number stream
Interpolation
1-D Interpolation

interp1	1-D data interpolation (table lookup)
griddedInterpolant	Gridded data interpolation
pchip	Piecewise Cubic Hermite Interpolating Polynomial (PCHIP)
spline	Cubic spline data interpolation
ppval	Evaluate piecewise polynomial
mkpp	Make piecewise polynomial
unmkpp	Piecewise polynomial details
padecoef	Padé approximation of time delays
interpft	1-D interpolation using FFT method
Gridded Data Interpolation

interp2	Interpolation for 2-D gridded data in meshgrid format
interp3	Interpolation for 3-D gridded data in meshgrid format
interpn	Interpolation for 1-D, 2-D, 3-D, and N-D gridded data in ndgrid format
griddedInterpolant	Gridded data interpolation
ndgrid	Rectangular grid in N-D space
meshgrid	Rectangular grid in 2-D and 3-D space
Scattered Data Interpolation

griddata	Interpolate scattered data
griddatan	Data gridding and hypersurface fitting (dimension ≥ 2)
scatteredInterpolant	Scattered data interpolation
Optimization
fminbnd	Find minimum of single-variable function on fixed interval
fminsearch	Find minimum of unconstrained multivariable function using derivative-free method
fzero	Root of nonlinear function
lsqnonneg	Solve nonnegative least-squares constraints problem
optimget	Optimization options values
optimset	Create or edit optimization options structure
Numerical Integration and Differential Equations
Ordinary Differential Equations

ode45	Solve nonstiff differential equations; medium order method
ode15s	Solve stiff differential equations and DAEs; variable order method
ode23	Solve nonstiff differential equations; low order method
ode113	Solve nonstiff differential equations; variable order method
ode23t	Solve moderately stiff ODEs and DAEs; trapezoidal rule
ode23tb	Solve stiff differential equations; low order method
ode23s	Solve stiff differential equations; low order method
ode15i	Solve fully implicit differential equations, variable order method
decic	Compute consistent initial conditions for ode15i
odextend	Extend solution of initial value problem for ordinary differential equation
odeget	Ordinary differential equation options parameters
odeset	Create or alter options structure for ordinary differential equation solvers
deval	Evaluate solution of differential equation problem
Boundary Value Problems

bvp4c	Solve boundary value problems for ordinary differential equations
bvp5c	Solve boundary value problems for ordinary differential equations
bvpinit	Form initial guess for BVP solvers
bvpxtend	Form guess structure for extending boundary value solutions
bvpget	Extract properties from options structure created with bvpset
bvpset	Create or alter options structure of boundary value problem
deval	Evaluate solution of differential equation problem
Delay Differential Equations

dde23	Solve delay differential equations (DDEs) with constant delays
ddesd	Solve delay differential equations (DDEs) with general delays
ddensd	Solve delay differential equations (DDEs) of neutral type
ddeget	Extract properties from delay differential equations options structure
ddeset	Create or alter delay differential equations options structure
deval	Evaluate solution of differential equation problem
Partial Differential Equations

pdepe	Solve initial-boundary value problems for parabolic-elliptic PDEs in 1-D
pdeval	Evaluate numerical solution of PDE using output of pdepe
Numerical Integration and Differentiation

integral	Numerically evaluate integral
integral2	Numerically evaluate double integral
integral3	Numerically evaluate triple integral
quadgk	Numerically evaluate integral, adaptive Gauss-Kronrod quadrature
quad2d	Numerically evaluate double integral, tiled method
cumtrapz	Cumulative trapezoidal numerical integration
trapz	Trapezoidal numerical integration
polyint	Integrate polynomial analytically
del2	Discrete Laplacian
diff	Differences and Approximate Derivatives
gradient	Numerical gradient
polyder	Polynomial derivative
Fourier Analysis and Filtering
abs	Absolute value and complex magnitude
angle	Phase angle
cplxpair	Sort complex numbers into complex conjugate pairs
fft	Fast Fourier transform
fft2	2-D fast Fourier transform
fftn	N-D fast Fourier transform
fftshift	Shift zero-frequency component to center of spectrum
fftw	Interface to FFTW library run-time algorithm tuning control
ifft	Inverse fast Fourier transform
ifft2	2-D inverse fast Fourier transform
ifftn	N-D inverse fast Fourier transform
ifftshift	Inverse FFT shift
nextpow2	Exponent of next higher power of 2
unwrap	Correct phase angles to produce smoother phase plots
conv	Convolution and polynomial multiplication
conv2	2-D convolution
convn	N-D convolution
deconv	Deconvolution and polynomial division
detrend	Remove linear trends
filter	1-D digital filter
filter2	2-D digital filter
Sparse Matrices
Sparse Matrix Creation

spdiags	Extract and create sparse band and diagonal matrices
speye	Sparse identity matrix
sprand	Sparse uniformly distributed random matrix
sprandn	Sparse normally distributed random matrix
sprandsym	Sparse symmetric random matrix
sparse	Create sparse matrix
spconvert	Import matrix from sparse matrix external format
Sparse Matrix Manipulation

issparse	Determine whether input is sparse
nnz	Number of nonzero matrix elements
nonzeros	Nonzero matrix elements
nzmax	Amount of storage allocated for nonzero matrix elements
spalloc	Allocate space for sparse matrix
spfun	Apply function to nonzero sparse matrix elements
spones	Replace nonzero sparse matrix elements with ones
spparms	Set parameters for sparse matrix routines
spy	Visualize sparsity pattern
find	Find indices and values of nonzero elements
full	Convert sparse matrix to full matrix
Reordering Algorithms

amd	Approximate minimum degree permutation
colamd	Column approximate minimum degree permutation
colperm	Sparse column permutation based on nonzero count
dmperm	Dulmage-Mendelsohn decomposition
randperm	Random permutation
symamd	Symmetric approximate minimum degree permutation
symrcm	Sparse reverse Cuthill-McKee ordering
Sparse Linear Algebra

condest	1-norm condition number estimate
eigs	Largest eigenvalues and eigenvectors of matrix
ichol	Incomplete Cholesky factorization
ilu	Sparse incomplete LU factorization
normest	2-norm estimate
spaugment	Form least squares augmented system
sprank	Structural rank
svds	Find singular values and vectors
Linear Equations (Iterative Methods)

bicg	Biconjugate gradients method
bicgstab	Biconjugate gradients stabilized method
bicgstabl	Biconjugate gradients stabilized (l) method
cgs	Conjugate gradients squared method
gmres	Generalized minimum residual method (with restarts)
lsqr	LSQR method
minres	Minimum residual method
pcg	Preconditioned conjugate gradients method
qmr	Quasi-minimal residual method
symmlq	Symmetric LQ method
tfqmr	Transpose-free quasi-minimal residual method
Graph and Tree Algorithms

etree	Elimination tree
etreeplot	Plot elimination tree
gplot	Plot nodes and links representing adjacency matrix
symbfact	Symbolic factorization analysis
treelayout	Lay out tree or forest
treeplot	Plot picture of tree
unmesh	Convert edge matrix to coordinate and Laplacian matrices
Computational Geometry
Triangulation Representation

triangulation	Triangulation in 2-D or 3-D
tetramesh	Tetrahedron mesh plot
trimesh	Triangular mesh plot
triplot	2-D triangular plot
trisurf	Triangular surface plot
Delaunay Triangulation

delaunayTriangulation	Delaunay triangulation in 2-D and 3-D
delaunay	Delaunay triangulation
delaunayn	N-D Delaunay triangulation
tetramesh	Tetrahedron mesh plot
trimesh	Triangular mesh plot
triplot	2-D triangular plot
trisurf	Triangular surface plot
Spatial Search

triangulation	Triangulation in 2-D or 3-D
delaunayTriangulation	Delaunay triangulation in 2-D and 3-D
dsearchn	N-D nearest point search
tsearchn	N-D closest simplex search
delaunay	Delaunay triangulation
delaunayn	N-D Delaunay triangulation
Convex Hull

convhull	Convex hull
convhulln	N-D convex hull
patch	Create one or more filled polygons
trisurf	Triangular surface plot
Voronoi Diagram

patch	Create one or more filled polygons
voronoi	Voronoi diagram
voronoin	N-D Voronoi diagram
Elementary Polygons

polyarea	Area of polygon
inpolygon	Points inside polygonal region
rectint	Rectangle intersection area
Graphics

2-D and 3-D Plots
Line Plots

plot	2-D line plot
plotyy	2-D line plots with y-axes on both left and right side
plot3	3-D line plot
loglog	Log-log scale plot
semilogx	Semilogarithmic plot
semilogy	Semilogarithmic plot
errorbar	Plot error bars along curve
fplot	Plot function between specified limits
ezplot	Easy-to-use function plotter
ezplot3	Easy-to-use 3-D parametric curve plotter
LineSpec (Line Specification)	Line specification string syntax
ColorSpec (Color Specification)	Color specification
Pie Charts, Bar Plots, and Histograms

bar	Bar graph
bar3	Plot 3-D bar graph
barh	Plot bar graph horizontally
bar3h	Plot horizontal 3-D bar graph
hist	Histogram plot
histc	Histogram bin count
rose	Angle histogram plot
pareto	Pareto chart
area	Filled area 2-D plot
pie	Pie chart
pie3	3-D pie chart
Discrete Data Plots

stem	Plot discrete sequence data
stairs	Stairstep graph
stem3	Plot 3-D discrete sequence data
scatter	Scatter plot
scatter3	3-D scatter plot
spy	Visualize sparsity pattern
plotmatrix	Scatter plot matrix
Polar Plots

polar	Polar coordinate plot
rose	Angle histogram plot
compass	Plot arrows emanating from origin
ezpolar	Easy-to-use polar coordinate plotter
LineSpec (Line Specification)	Line specification string syntax
ColorSpec (Color Specification)	Color specification
Contour Plots

contour	Contour plot of matrix
contourf	Filled 2-D contour plot
contourc	Low-level contour plot computation
contour3	3-D contour plot
contourslice	Draw contours in volume slice planes
ezcontour	Easy-to-use contour plotter
ezcontourf	Easy-to-use filled contour plotter
Vector Fields

feather	Plot velocity vectors
quiver	Quiver or velocity plot
compass	Plot arrows emanating from origin
quiver3	3-D quiver or velocity plot
streamslice	Plot streamlines in slice planes
streamline	Plot streamlines from 2-D or 3-D vector data
Surfaces, Volumes, and Polygons

Surface and Mesh Plots
surf	3-D shaded surface plot
surfc	Contour plot under a 3-D shaded surface plot
surface	Create surface object
surfl	Surface plot with colormap-based lighting
surfnorm	Compute and display 3-D surface normals
mesh	Mesh plot
meshc	Plot a contour graph under mesh graph
meshz	Plot a curtain around mesh plot
waterfall	Waterfall plot
ribbon	Ribbon plot
contour3	3-D contour plot
peaks	Example function of two variables
cylinder	Generate cylinder
ellipsoid	Generate ellipsoid
sphere	Generate sphere
pcolor	Pseudocolor (checkerboard) plot
surf2patch	Convert surface data to patch data
ezsurf	Easy-to-use 3-D colored surface plotter
ezsurfc	Easy-to-use combination surface/contour plotter
ezmesh	Easy-to-use 3-D mesh plotter
ezmeshc	Easy-to-use combination mesh/contour plotter
Volume Visualization
contourslice	Draw contours in volume slice planes
flow	Simple function of three variables
isocaps	Compute isosurface end-cap geometry
isocolors	Calculate isosurface and patch colors
isonormals	Compute normals of isosurface vertices
isosurface	Extract isosurface data from volume data
reducepatch	Reduce number of patch faces
reducevolume	Reduce number of elements in volume data set
shrinkfaces	Reduce size of patch faces
slice	Volumetric slice plot
smooth3	Smooth 3-D data
subvolume	Extract subset of volume data set
volumebounds	Coordinate and color limits for volume data
coneplot	Plot velocity vectors as cones in 3-D vector field
curl	Compute curl and angular velocity of vector field
divergence	Compute divergence of vector field
interpstreamspeed	Interpolate stream-line vertices from flow speed
stream2	Compute 2-D streamline data
stream3	Compute 3-D streamline data
streamline	Plot streamlines from 2-D or 3-D vector data
streamparticles	Plot stream particles
streamribbon	3-D stream ribbon plot from vector volume data
streamslice	Plot streamlines in slice planes
streamtube	Create 3-D stream tube plot
Polygons
fill	Filled 2-D polygons
fill3	Filled 3-D polygons
patch	Create one or more filled polygons
surf2patch	Convert surface data to patch data
Animation

movie	Play recorded movie frames
noanimate	Change EraseMode of all objects to normal
drawnow	Update figure window and execute pending callbacks
refreshdata	Refresh data in graph when data source is specified
frame2im	Return image data associated with movie frame
getframe	Capture movie frame
im2frame	Convert image to movie frame
comet	2-D comet plot
comet3	3-D comet plot
Formatting and Annotation
Titles and Labels

title	Add title to current axes
xlabel	Label x-axis
ylabel	Label y-axis
zlabel	Label z-axis
clabel	Contour plot elevation labels
datetick	Date formatted tick labels
texlabel	Format text into TeX string
legend	Graph legend for lines and patches
colorbar	Colorbar showing color scale
Coordinate System

xlim	Set or query x-axis limits
ylim	Set or query y-axis limits
zlim	Set or query z-axis limits
box	Axes border
grid	Grid lines for 2-D and 3-D plots
daspect	Set or query axes data aspect ratio
pbaspect	Set or query plot box aspect ratio
axes	Create axes graphics object
axis	Axis scaling and appearance
subplot	Create axes in tiled positions
hold	Retain current graph when adding new graphs
gca	Current axes handle
cla	Clear current axes
Annotation

annotation	Create annotation objects
text	Create text object in current axes
legend	Graph legend for lines and patches
title	Add title to current axes
xlabel	Label x-axis
ylabel	Label y-axis
zlabel	Label z-axis
datacursormode	Enable, disable, and manage interactive data cursor mode
ginput	Graphical input from mouse or cursor
gtext	Mouse placement of text in 2-D view
Colormaps

colormap	Set and get current colormap
colormapeditor	Open colormap editor
colorbar	Colorbar showing color scale
brighten	Brighten or darken colormap
contrast	Grayscale colormap for contrast enhancement
shading	Set color shading properties
graymon	Set default figure properties for grayscale monitors
caxis	Color axis scaling
hsv2rgb	Convert HSV colormap to RGB colormap
rgb2hsv	Convert RGB colormap to HSV colormap
rgbplot	Plot colormap
spinmap	Spin colormap
colordef	Set default property values to display different color schemes
whitebg	Change axes background color
Data Exploration

hidden	Remove hidden lines from mesh plot
pan	Pan view of graph interactively
reset	Reset graphics object properties to their defaults
rotate	Rotate object about specified origin and direction
rotate3d	Rotate 3-D view using mouse
selectmoveresize	Select, move, resize, or copy axes and uicontrol graphics objects
zoom	Turn zooming on or off or magnify by factorMagnify by a factor
datacursormode	Enable, disable, and manage interactive data cursor mode
figurepalette	Show or hide Figure Palette
plotbrowser	Show or hide figure Plot Browser
plotedit	Interactively edit and annotate plots
plottools	Show or hide plot tools
propertyeditor	Show or hide Property Editor
showplottool	Show or hide figure plot tool
Data Brushing

brush	Interactively mark, delete, modify, and save observations in graphs
datacursormode	Enable, disable, and manage interactive data cursor mode
linkdata	Automatically update graphs when variables change
refreshdata	Refresh data in graph when data source is specified
3-D Scene Control

Camera Views
view	Viewpoint specification
makehgtform	Create 4-by-4 transform matrix
viewmtx	View transformation matrices
cameratoolbar	Control camera toolbar programmatically
campan	Rotate camera target around camera position
camzoom	Zoom in and out on scene
camdolly	Move camera position and target
camlookat	Position camera to view object or group of objects
camorbit	Rotate camera position around camera target
campos	Set or query camera position
camproj	Set or query projection type
camroll	Rotate camera about view axis
camtarget	Set or query location of camera target
camup	Set or query camera up vector
camva	Set or query camera view angle
Lighting and Transparency
camlight	Create or move light object in camera coordinates
light	Create light object
lightangle	Create or position light object in spherical coordinates
lighting	Specify lighting algorithm
diffuse	Calculate diffuse reflectance
material	Control reflectance properties of surfaces and patches
specular	Calculate specular reflectance
alim	Set or query axes alpha limits
alpha	Set transparency properties for objects in current axes
alphamap	Specify figure alphamap (transparency)
Images
Image File Operations

image	Display image object
imagesc	Scale data and display image object
imread	Read image from graphics file
imwrite	Write image to graphics file
imfinfo	Information about graphics file
imformats	Manage image file format registry
frame2im	Return image data associated with movie frame
im2frame	Convert image to movie frame
im2java	Convert image to Java image
Modifying Images

ind2rgb	Convert indexed image to RGB image
rgb2ind	Convert RGB image to indexed image
imapprox	Approximate indexed image by reducing number of colors
dither	Convert image, increasing apparent color resolution by dithering
cmpermute	Rearrange colors in colormap
cmunique	Eliminate duplicate colors in colormap; convert grayscale or truecolor image to indexed image
Printing and Exporting
print	Print figure or save to file
printopt	Configure printer defaults
printdlg	Print dialog box
printpreview	Preview figure to print
orient	Hardcopy paper orientation
savefig	Save figure to FIG-file
openfig	Open new copy or raise existing copy of saved figure
hgexport	Export figure
hgsave	Save Handle Graphics object hierarchy to file
hgload	Load Handle Graphics object hierarchy from file
saveas	Save figure or Simulink block diagram using specified format
Graphics Objects
Graphics Object Identification

gca	Current axes handle
gcf	Current figure handle
gcbf	Handle of figure containing object whose callback is executing
gcbo	Handle of object whose callback is executing
gco	Handle of current object
ancestor	Ancestor of graphics object
allchild	Find all children of specified objects
findall	Find all graphics objects
findfigs	Find visible offscreen figures
findobj	Locate graphics objects with specific properties
gobjects	Create array of graphics handles
ishghandle	True for Handle Graphics object handles
ishandle	Test for valid graphics or Java object handle
copyobj	Copy graphics objects and their descendants
delete	Remove files or objects
get	Query Handle Graphics object properties
set	Set Handle Graphics object properties
propedit	Open Property Editor
Core Objects

root object	Root
figure	Create figure graphics object
axes	Create axes graphics object
image	Display image object
light	Create light object
line	Create line object
patch	Create one or more filled polygons
rectangle	Create 2-D rectangle object
surface	Create surface object
text	Create text object in current axes
Annotation Objects

annotation	Create annotation objects
Plot Objects

set	Set Handle Graphics object properties
get	Query Handle Graphics object properties
Group Objects

hggroup	Create hggroup object
hgtransform	Create hgtransform graphics object
makehgtform	Create 4-by-4 transform matrix
Figure Windows

figure	Create figure graphics object
gcf	Current figure handle
close	Remove specified figure
clf	Clear current figure window
refresh	Redraw current figure
newplot	Determine where to draw graphics objects
shg	Show most recent graph window
closereq	Default figure close request function
dragrect	Drag rectangles with mouse
drawnow	Update figure window and execute pending callbacks
rbbox	Create rubberband box for area selection
opengl	Control OpenGL rendering
Axes Property Operations

axes	Create axes graphics object
hold	Retain current graph when adding new graphs
ishold	Current hold state
newplot	Determine where to draw graphics objects
Object Property Operations

linkaxes	Synchronize limits of specified 2-D axes
linkprop	Keep same value for corresponding properties of graphics objects
refreshdata	Refresh data in graph when data source is specified
waitfor	Block execution and wait for event or condition
get	Query Handle Graphics object properties
set	Set Handle Graphics object properties
Programming Scripts and Functions

Control Flow
if, elseif, else	Execute statements if condition is true
for	Execute statements specified number of times
parfor	Parallel for loop
switch, case, otherwise	Switch among several cases based on expression
try, catch	Execute statements and catch resulting errors
while	Repeatedly execute statements while condition is true
break	Terminate execution of for or while loop
continue	Pass control to next iteration of for or while loop
end	Terminate block of code, or indicate last array index
pause	Halt execution temporarily
return	Return to invoking function
Scripts
edit	Edit or create file
input	Request user input
publish	Generate view of MATLAB file in specified format
notebook	Open MATLAB Notebook in Microsoft Word software (on Microsoft Windows platforms)
grabcode	Extract MATLAB code from file published to HTML
snapnow	Force snapshot of image for inclusion in published document
Functions
Function Basics

function	Declare function name, inputs, and outputs
Input and Output Arguments

nargin	Number of function input arguments
nargout	Number of function output arguments
varargin	Variable-length input argument list
varargout	Variable-length output argument list
narginchk	Validate number of input arguments
nargoutchk	Validate number of output arguments
validateattributes	Check validity of array
validatestring	Check validity of text string
inputParser	Parse function inputs
inputname	Variable name of function input
Variables

persistent	Define persistent variable
isvarname	Determine whether input is valid variable name
matlab.lang.makeUniqueStrings	Construct unique strings from input strings
matlab.lang.makeValidName	Construct valid MATLAB identifiers from input strings
namelengthmax	Maximum identifier length
assignin	Assign value to variable in specified workspace
global	Declare global variables
isglobal	Determine whether input is global variable
Error Handling

try, catch	Execute statements and catch resulting errors
error	Display message and abort function
warning	Warning message
lastwarn	Last warning message
assert	Generate error when condition is violated
onCleanup	Cleanup tasks upon function completion
Debugging
dbclear	Clear breakpoints
dbcont	Resume execution
dbdown	Reverse workspace shift performed by dbup, while in debug mode
dbquit	Quit debug mode
dbstack	Function call stack
dbstatus	List all breakpoints
dbstep	Execute one or more lines from current breakpoint
dbstop	Set breakpoints for debugging
dbtype	List text file with line numbers
dbup	Shift current workspace to workspace of caller, while in debug mode
checkcode	Check MATLAB code files for possible problems
keyboard	Input from keyboard
mlintrpt	Run checkcode for file or folder, reporting results in browser
Coding and Productivity Tips
edit	Edit or create file
Programming Utilities
echo	Display statements during function execution
eval	Execute MATLAB expression in text string
evalc	Evaluate MATLAB expression with capture
evalin	Execute MATLAB expression in specified workspace
feval	Evaluate function
run	Run MATLAB script
builtin	Execute built-in function from overloaded method
matlab.codetools.requiredFilesAndProducts	List dependencies of MATLAB program files
mfilename	File name of currently running function
pcode	Create protected function file
timer	Create object to schedule execution of MATLAB commands
Data and File Management

Workspace Variables
clear	Remove items from workspace, freeing up system memory
clearvars	Clear variables from memory
disp	Display text or array
openvar	Open workspace variable in Variables editor or other graphical editing tool
who	List variables in workspace
whos	List variables in workspace, with sizes and types
load	Load variables from file into workspace
save	Save workspace variables to file
matfile	Access and change variables directly in MAT-files, without loading into memory
Data Import and Export
Import and Export Basics

importdata	Load data from file
uiimport	Import data interactively
Text Files

csvread	Read comma-separated value file
csvwrite	Write comma-separated value file
dlmread	Read ASCII-delimited file of numeric data into matrix
dlmwrite	Write matrix to ASCII-delimited file
textscan	Read formatted data from text file or string
readtable	Create table from file
writetable	Write table to file
type	Display contents of file
Spreadsheets

xlsfinfo	Determine if file contains Microsoft Excel spreadsheet
xlsread	Read Microsoft Excel spreadsheet file
xlswrite	Write Microsoft Excel spreadsheet file
readtable	Create table from file
writetable	Write table to file
Low-Level File I/O

fclose	Close one or all open files
feof	Test for end-of-file
ferror	Information about file I/O errors
fgetl	Read line from file, removing newline characters
fgets	Read line from file, keeping newline characters
fileread	Read contents of file into string
fopen	Open file, or obtain information about open files
fprintf	Write data to text file
fread	Read data from binary file
frewind	Move file position indicator to beginning of open file
fscanf	Read data from text file
fseek	Move to specified position in file
ftell	Position in open file
fwrite	Write data to binary file
Images

im2java	Convert image to Java image
imfinfo	Information about graphics file
imread	Read image from graphics file
imwrite	Write image to graphics file
Tiff	MATLAB Gateway to LibTIFF library routines
Scientific Data

netCDF Files
nccreate	Create variable in NetCDF file
ncdisp	Display contents of NetCDF data source in Command Window
ncinfo	Return information about NetCDF data source
ncread	Read data from variable in NetCDF data source
ncreadatt	Read attribute value from NetCDF data source
ncwrite	Write data to NetCDF file
ncwriteatt	Write attribute to NetCDF file
ncwriteschema	Add NetCDF schema definitions to NetCDF file
netcdf	Summary of MATLAB Network Common Data Form (NetCDF) capabilities
HDF5 Files
High-Level Functions
h5create	Create HDF5 data set
h5disp	Display contents of HDF5 file
h5info	Return information about HDF5 file
h5read	Read data from HDF5 data set
h5readatt	Read attribute from HDF5 file
h5write	Write to HDF5 data set
h5writeatt	Write HDF5 attribute
Low-Level Functions
Library (H5)
H5.close	Close HDF5 library
H5.garbage_collect	Free unused memory in HDF5 library
H5.get_libversion	Version of HDF5 library
H5.open	Open HDF5 library
H5.set_free_list_limits	Set size limits on free lists
Attribute (H5A)
H5A.close	Close specified attribute
H5A.create	Create attribute
H5A.delete	Delete attribute
H5A.get_info	Information about attribute
H5A.get_name	Attribute name
H5A.get_space	Copy of attribute data space
H5A.get_type	Copy of attribute data type
H5A.iterate	Execute function for attributes attached to object
H5A.open	Open attribute
H5A.open_by_idx	Open attribute specified by index
H5A.open_by_name	Open attribute specified by name
H5A.read	Read attribute
H5A.write	Write attribute
Dataset (H5D)
H5D.close	Close dataset
H5D.create	Create new dataset
H5D.get_access_plist	Copy of dataset access property list
H5D.get_create_plist	Copy of dataset creation property list
H5D.get_offset	Location of dataset in file
H5D.get_space	Copy of dataset data space
H5D.get_space_status	Determine if space is allocated
H5D.get_storage_size	Determine required storage size
H5D.get_type	Copy of datatype
H5D.open	Open specified dataset
H5D.read	Read data from HDF5 dataset
H5D.set_extent	Change size of dataset dimensions
H5D.vlen_get_buf_size	Determine variable length storage requirements
H5D.write	Write data to HDF5 dataset
Dimension Scale (H5DS)
H5DS.attach_scale	Attach dimension scale to specific dataset dimension
H5DS.detach_scale	Detach dimension scale from specific dataset dimension
H5DS.get_label	Retrieve label from specific dataset dimension
H5DS.get_num_scales	Number of scales attached to dataset dimension
H5DS.get_scale_name	Name of dimension scale
H5DS.is_scale	Determine if dataset is a dimension scale
H5DS.iterate_scales	Iterate on scales attached to dataset dimension
H5DS.set_label	Set label for dataset dimension
H5DS.set_scale	Convert dataset to dimension scale
Error (H5E)
H5E.clear	Clear error stack
H5E.get_major	Description of major error number
H5E.get_minor	Description of minor error number
H5E.walk	Walk error stack
File (H5F)
H5F.close	Close HDF5 file
H5F.create	Create HDF5 file
H5F.flush	Flush buffers to disk
H5F.get_access_plist	File access property list
H5F.get_create_plist	File creation property list
H5F.get_filesize	Size of HDF5 file
H5F.get_freespace	Amount of free space in file
H5F.get_info	Global information about file
H5F.get_mdc_config	Metadata cache configuration
H5F.get_mdc_hit_rate	Metadata cache hit-rate
H5F.get_mdc_size	Metadata cache size data
H5F.get_name	Name of HDF5 file
H5F.get_obj_count	Number of open objects in HDF5 file
H5F.get_obj_ids	List of open HDF5 file objects
H5F.is_hdf5	Determine if file is HDF5
H5F.mount	Mount HDF5 file onto specified location
H5F.open	Open HDF5 file
H5F.reopen	Reopen HDF5 file
H5F.set_mdc_config	Configure HDF5 file metadata cache
H5F.unmount	Unmount file or group from mount point
Group (H5G)
H5G.close	Close group
H5G.create	Create group
H5G.get_info	Information about group
H5G.open	Open specified group
Identifier (H5I)
H5I.dec_ref	Decrement reference count
H5I.get_file_id	File identifier for specified object
H5I.get_name	Name of object
H5I.get_ref	Reference count of object
H5I.get_type	Type of object
H5I.inc_ref	Increment reference count of specified object
H5I.is_valid	Determine if specified identifier is valid
Link (H5L)
H5L.copy	Copy link from source location to destination location
H5L.create_external	Create soft link to external object
H5L.create_hard	Create hard link
H5L.create_soft	Create soft link
H5L.delete	Remove link
H5L.exists	Determine if link exists
H5L.get_info	Information about link
H5L.get_name_by_idx	Information about link specified by index
H5L.get_val	Value of symbolic link
H5L.iterate	Iterate over links
H5L.iterate_by_name	Iterate through links in group specified by name
H5L.move	Rename link
H5L.visit	Recursively iterate through links in group specified by group identifier
H5L.visit_by_name	Recursively iterate through links in group specified by location and group name
MATLAB (H5ML)
H5ML.compare_values	Numerically compare two HDF5 values
H5ML.get_constant_names	Constants known by HDF5 library
H5ML.get_constant_value	Value corresponding to a string
H5ML.get_function_names	Functions provided by HDF5 library
H5ML.get_mem_datatype	Data type for dataset ID
Object (H5O)
H5O.close	Close object
H5O.copy	Copy object from source location to destination location
H5O.get_comment	Get comment for object specified by object identifier
H5O.get_comment_by_name	Get comment for object specified by location and object name
H5O.get_info	Object metadata
H5O.link	Create hard link to specified object
H5O.open	Open specified object
H5O.open_by_idx	Open object specified by index
H5O.set_comment	Set comment for object specified by object identifier
H5O.set_comment_by_name	Set comment for object specified by location and object name
H5O.visit	Visit objects specified by object identifier
H5O.visit_by_name	Visit objects specified by location and object name
Property (H5P) General Property List Operations
H5P.close	Close property list
H5P.copy	Copy of property list
H5P.create	Create new property list
H5P.get_class	Property list class
Generic Property List Operations
H5P.close_class	Close property list class
H5P.equal	Determine equality of property lists
H5P.exist	Determine if specified property exists in property list
H5P.get	Value of specified property in property list
H5P.get_class_name	Name of property list class
H5P.get_class_parent	Identifier for parent class
H5P.get_nprops	Query number of properties in property list or class
H5P.get_size	Query size of property value in bytes
H5P.isa_class	Determine if property list is member of class
H5P.iterate	Iterate over properties in property list
H5P.set	Set property list value
Dataset Access, Memory, and Transfer Properties
H5P.get_btree_ratios	B-tree split ratios
H5P.get_chunk_cache	Raw data chunk cache parameters
H5P.get_dxpl_multi	Data access property lists for multiple files
H5P.get_edc_check	Determine if error detection is enabled
H5P.get_hyper_vector_size	Number of I/O vectors
H5P.set_btree_ratios	Set B-tree split ratios for dataset transfer
H5P.set_chunk_cache	Set raw data chunk cache parameters
H5P.set_dxpl_multi	Set data transfer property list for multifile driver
H5P.set_edc_check	Enable error detection for dataset transfer
H5P.set_hyper_vector_size	Set number of I/O vectors for hyperslab I/O
Dataset Creation Properties
H5P.all_filters_avail	Determine availability of all filters
H5P.fill_value_defined	Determine if fill value is defined
H5P.get_alloc_time	Return timing of storage space allocation
H5P.get_chunk	Return size of chunks
H5P.get_external	Return information about external file
H5P.get_external_count	Return count of external files
H5P.get_fill_time	Return time when fill values are written to dataset
H5P.get_fill_value	Return dataset fill value
H5P.get_filter	Return information about filter in pipeline
H5P.get_filter_by_id	Return information about specified filter
H5P.get_layout	Determine layout of raw data for dataset
H5P.get_nfilters	Return number of filters in pipeline
H5P.modify_filter	Modify filter in pipeline
H5P.remove_filter	Remove filter from property list
H5P.set_alloc_time	Set timing for storage space allocation
H5P.set_chunk	Set chunk size
H5P.set_deflate	Set compression method and compression level
H5P.set_external	Add additional file to external file list
H5P.set_fill_time	Set time when fill values are written to dataset
H5P.set_fill_value	Set fill value for dataset creation property list
H5P.set_filter	Add filter to filter pipeline
H5P.set_fletcher32	Set Fletcher32 checksum filter in dataset creation
H5P.set_layout	Set type of storage for dataset
H5P.set_nbit	Set N-Bit filter
H5P.set_scaleoffset	Set Scale-Offset filter
H5P.set_shuffle	Set shuffle filter
File Access Properties
H5P.get_alignment	Retrieve alignment properties
H5P.get_driver	Low-level file driver
H5P.get_family_offset	Offset for family file driver
H5P.get_fapl_core	Information about core file driver properties
H5P.get_fapl_family	File access property list information
H5P.get_fapl_multi	Information about multifile access property list
H5P.get_fclose_degree	File close degree
H5P.get_libver_bounds	Library version bounds settings
H5P.get_gc_references	Garbage collection references setting
H5P.get_mdc_config	Metadata cache configuration
H5P.get_meta_block_size	Metadata block size setting
H5P.get_multi_type	Type of data property for MULTI driver
H5P.get_sieve_buf_size	Maximum data sieve buffer size
H5P.get_small_data_block_size	Small data block size setting
H5P.set_alignment	Set alignment properties for file access property list
H5P.set_family_offset	Set offset property for family of files
H5P.set_fapl_core	Modify file access to use H5FD_CORE driver
H5P.set_fapl_family	Set file access to use family driver
H5P.set_fapl_log	Set use of logging driver
H5P.set_fapl_multi	Set use of multifile driver
H5P.set_fapl_sec2	Set file access for sec2 driver
H5P.set_fapl_split	Set file access for emulation of split file driver
H5P.set_fapl_stdio	Set file access for standard I/O driver
H5P.set_fclose_degree	Set file access for file close degree
H5P.set_gc_references	Set garbage collection references flag
H5P.set_libver_bounds	Set library version bounds for objects
H5P.set_mdc_config	Set initial metadata cache configuration
H5P.set_meta_block_size	Set minimum metadata block size
H5P.set_multi_type	Specify type of data accessed with MULTI driver
H5P.set_sieve_buf_size	Set maximum size of data sieve buffer
H5P.set_small_data_block_size	Set size of block reserved for small data
File Creation Properties
H5P.get_istore_k	Return 1/2 rank of indexed storage B-tree
H5P.get_sizes	Return size of offsets and lengths
H5P.get_sym_k	Return size of B-tree 1/2 rank and leaf node 1/2 size
H5P.get_userblock	Return size of user block
H5P.get_version	Return version information for file creation property list
H5P.set_istore_k	Set size of parameter for indexing chunked datasets
H5P.set_sizes	Set byte size of offsets and lengths
H5P.set_sym_k	Set size of parameters used to control symbol table nodes
H5P.set_userblock	Set user block size
Object Copy and Object Creation Properties
H5P.get_attr_creation_order	Return tracking order and indexing settings
H5P.get_attr_phase_change	Retrieve attribute phase change thresholds
H5P.get_copy_object	Return properties to be used when object is copied
H5P.set_attr_creation_order	Set tracking of attribute creation order
H5P.set_attr_phase_change	Set attribute storage phase change thresholds
H5P.set_copy_object	Set properties to be used when objects are copied
Group Creation Properties
H5P.get_create_intermediate_group	Determine creation of intermediate groups
H5P.get_link_creation_order	Query if link creation order is tracked
H5P.get_link_phase_change	Query settings for conversion between groups
H5P.set_create_intermediate_group	Set creation of intermediate groups
H5P.set_link_creation_order	Set creation order tracking and indexing
H5P.set_link_phase_change	Set parameters for group conversion
HDF5 String Properties
H5P.get_char_encoding	Return character encoding
H5P.set_char_encoding	Set character encoding used to encode strings
Reference (H5R)
H5R.create	Create reference
H5R.dereference	Open object specified by reference
H5R.get_name	Name of referenced object
H5R.get_obj_type	Type of referenced object
H5R.get_region	Copy of data space of specified region
Dataspace (H5S)
H5S.copy	Create copy of data space
H5S.create	Create new data space
H5S.close	Close data space
H5S.create_simple	Create new simple data space
H5S.extent_copy	Copy extent from source to destination data space
H5S.is_simple	Determine if data space is simple
H5S.offset_simple	Set offset of simple data space
H5S.select_all	Select entire extent of data space
H5S.select_elements	Specify coordinates to include in selection
H5S.select_hyperslab	Select hyperslab region
H5S.select_none	Reset selection region to include no elements
H5S.select_valid	Determine validity of selection
H5S.set_extent_none	Remove extent from data space
H5S.set_extent_simple	Set size of data space
H5S.get_select_bounds	Bounding box of data space selection
H5S.get_select_elem_npoints	Number of element points in selection
H5S.get_select_elem_pointlist	Element points in data space selection
H5S.get_select_hyper_blocklist	List of hyperslab blocks
H5S.get_select_hyper_nblocks	Number of hyperslab blocks
H5S.get_select_npoints	Number of elements in data space selection
H5S.get_select_type	Type of data space selection
H5S.get_simple_extent_dims	Data space size and maximum size
H5S.get_simple_extent_ndims	Data space rank
H5S.get_simple_extent_npoints	Number of elements in data space
H5S.get_simple_extent_type	Data space class
Datatype (H5T) General Data Type Operation
H5T.close	Close data type
H5T.commit	Commit transient data type
H5T.committed	Determine if data type is committed
H5T.copy	Copy data type
H5T.create	Create new data type
H5T.detect_class	Determine of data type contains specific class
H5T.equal	Determine equality of data types
H5T.get_class	Data type class identifier
H5T.get_create_plist	Copy of data type creation property list
H5T.get_native_type	Native data type of dataset data type
H5T.get_size	Size of data type in bytes
H5T.get_super	Base data type
H5T.lock	Lock data type
H5T.open	Open named data type
Array Data Type
H5T.array_create	Create array data type object
H5T.get_array_dims	Array dimension extents
H5T.get_array_ndims	Rank of array data type
Atomic Data Type Properties
H5T.get_cset	Character set of string data type
H5T.get_ebias	Exponent bias of floating-point type
H5T.get_fields	Floating-point data type bit field information
H5T.get_inpad	Internal padding type for floating-point data types
H5T.get_norm	Mantissa normalization type
H5T.get_offset	Bit offset of first significant bit
H5T.get_order	Byte order of atomic data type
H5T.get_pad	Padding type of least and most-significant bits
H5T.get_precision	Precision of atomic data type
H5T.get_sign	Sign type for integer data type
H5T.get_strpad	Storage mechanism for string data type
H5T.set_cset	Set character dataset for string data type
H5T.set_ebias	Set exponent bias of floating-point data type
H5T.set_fields	Set sizes and locations of floating-point bit fields
H5T.set_inpad	Specify how unused internal bits are to be filled
H5T.set_norm	Set mantissa normalization of floating-point data type
H5T.set_offset	Set bit offset of first significant bit
H5T.set_order	Set byte ordering of atomic data type
H5T.set_pad	Set padding type for least and most significant bits
H5T.set_precision	Set precision of atomic data type
H5T.set_sign	Set sign property for integer data type
H5T.set_size	Set size of data type in bytes
H5T.set_strpad	Set storage mechanism for string data type
Compound Data Type
H5T.get_member_class	Data type class for compound data type member
H5T.get_member_index	Index of compound or enumeration type member
H5T.get_member_name	Name of compound or enumeration type member
H5T.get_member_offset	Offset of field of compound data type
H5T.get_member_type	Data type of specified member
H5T.get_nmembers	Number of elements in enumeration type
H5T.insert	Add member to compound data type
H5T.pack	Recursively remove padding from compound data type
Enumeration Data Type
H5T.enum_create	Create new enumeration data type
H5T.enum_insert	Insert enumeration data type member
H5T.enum_nameof	Name of enumeration data type member
H5T.enum_valueof	Value of enumeration data type member
H5T.get_member_value	Value of enumeration data type member
Opaque Data Type Properties
H5T.get_tag	Tag associated with opaque data type
H5T.set_tag	Tag opaque data type with description
Variable-length Data Type
H5T.is_variable_str	Determine if data type is variable-length string
H5T.vlen_create	Create new variable-length data type
Filters and Compression (H5Z)
H5Z.filter_avail	Determine if filter is available
H5Z.get_filter_info	Information about filter
HDF4 Files
High-Level Functions
hdfinfo	Information about HDF4 or HDF-EOS file
hdfread	Read data from HDF4 or HDF-EOS file
imread	Read image from graphics file
imwrite	Write image to graphics file
Low-Level Functions
hdfan	Gateway to HDF multifile annotation (AN) interface
hdfhx	Gateway to HDF external data (HX) interface
hdfh	Gateway to HDF H interface
hdfhd	Gateway to HDF HD interface
hdfhe	Gateway to HDF HE interface
hdfml	Utilities for working with MATLAB HDF gateway functions
hdfpt	Interface to HDF-EOS Point object
hdfv	Gateway to HDF Vgroup (V) interface
hdfvf	Gateway to VF functions in HDF Vdata interface
hdfvh	Gateway to VH functions in HDF Vdata interface
hdfvs	Gateway to VS functions in HDF Vdata interface
hdfdf24	Gateway to HDF 24-bit raster image (DF24) interface
hdfdfr8	Gateway to HDF 8-bit raster image (DFR8) interface
FITS Files
High-Level Functions
fitsdisp	Display FITS metadata
fitsinfo	Information about FITS file
fitsread	Read data from FITS file
fitswrite	Write image to FITS file
Low-Level Functions
File Access
createFile	Create FITS file
openFile	Open FITS file
closeFile	Close FITS file
deleteFile	Delete FITS file
fileName	Name of FITS file
fileMode	I/O mode of FITS file
Image Manipulation
createImg	Create FITS image
getImgSize	Size of image
getImgType	Data type of image
insertImg	Insert FITS image after current image
readImg	Read image data
setBscale	Reset image scaling
writeImg	Write to FITS image
Keywords
readCard	Header record of keyword
readKey	Keyword
readKeyCmplx	Keyword as complex scalar value
readKeyDbl	Keyword as double precision value
readKeyLongLong	Keyword as int64
readKeyLongStr	Long string value
readKeyUnit	Physical units string from keyword
readRecord	Header record specified by number
writeComment	Write or append COMMENT keyword to CHU
writeDate	Write DATE keyword to CHU
writeKey	Update or add new keyword into current HDU
writeKeyUnit	Write physical units string
writeHistory	Write or append HISTORY keyword to CHU
deleteKey	Delete key by name
deleteRecord	Delete key by record number
getHdrSpace	Number of keywords in header
Header Data Unit (HDU) Access
copyHDU	Copy current HDU from one file to another
getHDUnum	Number of current HDU in FITS file
getHDUtype	Type of current HDU
getNumHDUs	Total number of HDUs in FITS file
movAbsHDU	Move to absolute HDU number
movNamHDU	Move to first HDU having specific type and keyword values
movRelHDU	Move relative number of HDUs from current HDU
writeChecksum	Compute and write checksum for current HDU
deleteHDU	Delete current HDU in FITS file
Image Compression
imgCompress	Compress HDU from one file into another
isCompressedImg	Determine if current image is compressed
setCompressionType	Set image compression type
setHCompScale	Set scale parameter for HCOMPRESS algorithm
setHCompSmooth	Set smoothing for images compressed with HCOMPRESS
setTileDim	Set tile dimensions
ASCII and Binary Tables
createTbl	Create new ASCII or binary table extension
insertCol	Insert column into table
insertRows	Insert rows into table
insertATbl	Insert ASCII table after current HDU
insertBTbl	Insert binary table after current HDU
deleteCol	Delete column from table
deleteRows	Delete rows from table
getAColParms	ASCII table information
getBColParms	Binary table information
getColName	Table column name
getColType	Scaled column data type, repeat value, width
getEqColType	Column data type, repeat value, width
getNumCols	Number of columns in table
getNumRows	Number of rows in table
readATblHdr	Read header information from current ASCII table
readBTblHdr	Read header information from current binary table
readCol	Read rows of ASCII or binary table column
setTscale	Reset image scaling
writeCol	Write elements into ASCII or binary table column
Utilities
getConstantValue	Numeric value of named constant
getVersion	Revision number of the CFITSIO library
getOpenFiles	List of open FITS files
Band-Interleaved Files
multibandread	Read band-interleaved data from binary file
multibandwrite	Write band-interleaved data to file
Common Data Format
cdfepoch	Convert MATLAB formatted dates to CDF formatted dates
cdfinfo	Information about Common Data Format (CDF) file
cdfread	Read data from Common Data Format (CDF) file
cdfwrite	Write data to Common Data Format (CDF) file
todatenum	Convert CDF epoch object to MATLAB datenum
cdflib	Summary of Common Data Format (CDF) capabilities
Audio and Video

Reading and Writing Files
audioinfo	Information about audio file
audioread	Read audio file
audiowrite	Write audio file
mmfileinfo	Information about multimedia file
movie2avi	Create Audio/Video Interleaved (AVI) file from MATLAB movie
VideoReader	Read video files
VideoWriter	Write video files
Recording and Playback
audiodevinfo	Information about audio device
audioplayer	Create object for playing audio
audiorecorder	Create object for recording audio
sound	Convert matrix of signal data to sound
soundsc	Scale data and play as sound
Utilities
beep	Produce operating system beep sound
lin2mu	Convert linear audio signal to mu-law
mu2lin	Convert mu-law audio signal to linear
XML Documents

xmlread	Read XML document and return Document Object Model node
xmlwrite	Write XML Document Object Model node
xslt	Transform XML document using XSLT engine
Memory Mapping

memmapfile	Create memory map to a file
File Operations
Files and Folders

dir	List folder contents
ls	List folder contents
pwd	Identify current folder
fileattrib	Set or get attributes of file or folder
exist	Check existence of variable, function, folder, or class
isdir	Determine whether input is folder
type	Display contents of file
visdiff	Compare two text files, MAT-Files, binary files, Zip files, or folders
what	List MATLAB files in folder
which	Locate functions and files
cd	Change current folder
copyfile	Copy file or folder
delete	Remove files or objects
recycle	Set option to move deleted files to recycle folder
mkdir	Make new folder
movefile	Move file or folder
rmdir	Remove folder
open	Open file in appropriate application
winopen	Open file in appropriate application (Windows)
File Name Construction

fileparts	Parts of file name and path
fullfile	Build full file name from parts
filemarker	Character to separate file name and internal function name
filesep	File separator for current platform
tempdir	Name of system's temporary folder
tempname	Unique name for temporary file
matlabroot	Root folder
toolboxdir	Root folder for specified toolbox
File Compression

zip	Compress files into zip file
unzip	Extract contents of zip file
gzip	Compress files into GNU zip files
gunzip	Uncompress GNU zip files
tar	Compress files into tar file
untar	Extract contents of tar file
Search Path
addpath	Add folders to search path
rmpath	Remove folders from search path
path	View or change search path
savepath	Save current search path
userpath	View or change user portion of search path
genpath	Generate path string
pathsep	Search path separator for current platform
pathtool	Open Set Path dialog box to view and change search path
restoredefaultpath	Restore default search path
Operating System Commands
clipboard	Copy and paste strings to and from system clipboard
computer	Information about computer on which MATLAB software is running
dos	Execute DOS command and return output
getenv	Environment variable
perl	Call Perl script using appropriate operating system executable
setenv	Set environment variable
system	Execute operating system command and return output
unix	Execute UNIX command and return output
winqueryreg	Item from Windows registry
Internet File Access
ftp	Connect to FTP server
sendmail	Send email message to address list
urlread	Download URL content to MATLAB string
urlwrite	Download URL content and save to file
web	Open Web page or file in browser
Serial Port Devices
delete (serial)	Remove serial port object from memory
fclose (serial)	Disconnect serial port object from device
fgetl (serial)	Read line of text from device and discard terminator
fgets (serial)	Read line of text from device and include terminator
fopen (serial)	Connect serial port object to device
fprintf (serial)	Write text to device
fread (serial)	Read binary data from device
fscanf (serial)	Read data from device, and format as text
fwrite (serial)	Write binary data to device
get (serial)	Serial port object properties
instrcallback	Event information when event occurs
instrfind	Read serial port objects from memory to MATLAB workspace
instrfindall	Find visible and hidden serial port objects
isvalid (serial)	Determine whether serial port objects are valid
readasync	Read data asynchronously from device
record	Record data and event information to file
serial	Create serial port object
serialbreak	Send break to device connected to serial port
set (serial)	Configure or display serial port object properties
stopasync	Stop asynchronous read and write operations
clear (serial)	Remove serial port object from MATLAB workspace
load (serial)	Load serial port objects and variables into MATLAB workspace
save (serial)	Save serial port objects and variables to file
disp (serial)	Serial port object summary information
length (serial)	Length of serial port object array
size (serial)	Size of serial port object array
Hardware Support
Raspberry Pi Hardware

Install and Set Up Support for Raspberry Pi Hardware
supportPackageInstaller	Find and install support for third-party hardware or software
raspi_examples	Open featured examples for this support package
targetupdater	Open Support Package Installer and update firmware on third-party hardware
matlabshared.supportpkg.checkForUpdate	Get information about installed support packages
matlabshared.supportpkg.getInstalled	Get information about installed support packages
Create a Connection to Raspberry Pi Hardware
raspi	Create connection to Raspberry Pi hardware
LEDs
raspi	Create connection to Raspberry Pi hardware
writeLED	Turn LED on or off
showLEDs	Show location, name, and color of user-controllable LEDs
GPIO Pins
raspi	Create connection to Raspberry Pi hardware
configureDigitalPin	Configure GPIO pin as input or output
readDigitalPin	Read logical value from GPIO input pin
writeDigitalPin	Write logical value to GPIO output pin
showPins	Show diagram of GPIO pins
Serial Port
raspi	Create connection to Raspberry Pi hardware
serialdev	Create connection to serial device
read	Read data from serial device
write	Write data to serial device
I2C Interface
raspi	Create connection to Raspberry Pi hardware
scanI2CBus	Scan I2C bus device addresses
i2cdev	Create connection to I2C device
read	Read data from I2C device
write	Write data to I2C device
readRegister	Read from register on I2C device
writeRegister	Write to register on I2C device
enableI2C	Enable I2C interface
disableI2C	Disable I2C interface
SPI Interface
raspi	Create connection to Raspberry Pi hardware
spidev	Create connection to SPI device
writeRead	Write data to and read data from SPI device
enableSPI	Enable SPI interface
disableSPI	Disable SPI interface
Camera Board
raspi	Create connection to Raspberry Pi hardware
cameraboard	Create connection to Raspberry Pi Camera Board Module
snapshot	Capture RGB image from Camera Board
record	Record video from Camera Board
stop	Stop video recording from Camera Board
Linux
raspi	Create connection to Raspberry Pi hardware
system	Run command in Linux shell on Raspberry Pi hardware
openShell	Open terminal on host computer to use Linux shell on Raspberry Pi hardware
getFile	Transfer file from Raspberry Pi hardware to host computer
putFile	Transfer file from host computer to target hardware
deleteFile	Delete file on target hardware
Webcams

webcamlist	List of Webcams connected to your system
webcam	Create webcam object to acquire images from Webcams
preview	Preview live video data from Webcam
snapshot	Acquire single image frame from a Webcam
closepreview	Close Webcam preview window
GUI Building

GUI Building Basics
guide	Open GUI Layout Editor
inspect	Open Property Inspector
Component Selection
GUI Controls and Indicators

figure	Create figure graphics object
axes	Create axes graphics object
uicontrol	Create user interface control object
uitable	Create 2-D graphic table GUI component
uipanel	Create panel container object
uibuttongroup	Create container object to exclusively manage radio buttons and toggle buttons
actxcontrol	Create Microsoft ActiveX control in figure window
Menus and Toolbars

uimenu	Create menus and menu items on figure windows
uicontextmenu	Create context menu
uitoolbar	Create toolbar on figure
uipushtool	Create push button on toolbar
uitoggletool	Create toggle button on toolbar
Predefined Dialog Boxes

dialog	Create and display empty dialog box
errordlg	Create and open error dialog box
helpdlg	Create and open help dialog box
msgbox	Create and open message dialog box
questdlg	Create and open question dialog box
uigetpref	Specify and conditionally open dialog box according to user preference
uisetpref	Manage preferences used in uigetpref
waitbar	Open or update wait bar dialog box
warndlg	Open warning dialog box
export2wsdlg	Export variables to workspace
inputdlg	Create and open input dialog box
listdlg	Create and open list-selection dialog box
uisetcolor	Open standard dialog box for setting object's color specification (ColorSpec)
uisetfont	Open standard dialog box for setting object's font characteristics
printdlg	Print dialog box
printpreview	Preview figure to print
uigetdir	Open standard dialog box for selecting directory
uigetfile	Open standard dialog box for retrieving files
uiopen	Interactively select file to open and load data
uiputfile	Open standard dialog box for saving files
uisave	Interactively save workspace variables to MAT-file
menu	Generate menu of choices for user input
Component Layout
align	Align user interface controls (uicontrols) and axes
movegui	Move GUI figure to specified location on screen
getpixelposition	Get component position in pixels
setpixelposition	Set component position in pixels
listfonts	List available system fonts
textwrap	Wrapped string matrix for given uicontrol
uistack	Reorder visual stacking order of objects
Coding GUI Behavior
uiwait	Block program execution and wait to resume
uiresume	Resume execution of blocked program
waitfor	Block execution and wait for event or condition
waitforbuttonpress	Wait for key press or mouse-button click
getappdata	Value of application-defined data
setappdata	Specify application-defined data
isappdata	True if application-defined data exists
rmappdata	Remove application-defined data
guidata	Store or retrieve GUI data
guihandles	Create structure of handles
Packaging GUIs as Apps
matlab.apputil.create	Create or modify app project file for packaging app into .mlappinstall file using interactive dialog box
matlab.apputil.package	Package app files into .mlappinstall file
matlab.apputil.install	Install app from a .mlappinstall file
matlab.apputil.run	Run app programmatically
matlab.apputil.getInstalledAppInfo	List installed app information
matlab.apputil.uninstall	Uninstall app
Advanced Software Development

Object-Oriented Programming
Class Syntax Fundamentals

classdef	Class definition keywords
class	Determine class of object
isa	Determine if input is object of specified class
isequal	Determine array equality
isobject	Determine if input is MATLAB object
enumeration	Display class enumeration members and names
events	Event names
methods	Class method names
properties	Class property names
Defining MATLAB Classes

Class Definition and Organization
classdef	Class definition keywords
import	Add package or class to current import list
Properties
properties	Class property names
isprop	Determine if property of object
dynamicprops	Abstract class used to derive handle class with dynamic properties
meta.property	meta.property class describes MATLAB class properties
Methods
methods	Class method names
ismethod	Determine if method of object
meta.method	meta.method class describes MATLAB class methods
Handle Classes
handle	Abstract class for deriving handle classes
hgsetget	Abstract class used to derive handle class with set and get methods
dynamicprops	Abstract class used to derive handle class with dynamic properties
matlab.mixin.Copyable	Superclass providing copy functionality for handle objects
delete	Handle object destructor
findobj	Find handle objects matching specified conditions
isa	Determine if input is object of specified class
isvalid	Is object valid handle class object
findprop	Find meta.property object associated with property name
relationaloperators	Equality and sorting of handle objects
Events
events	Event names
notify	Notify listeners that event is occurring
addlistener	Create event listener
event.EventData	Base class for all data objects passed to event listeners
event.listener	Class defining listener objects
event.PropertyEvent	Data for property events
event.proplistener	Define listener object for property events
Object Arrays
empty	Create empty array
matlab.mixin.Heterogeneous	Superclass for heterogeneous array formation
Class Hierarchies
superclasses	Superclass names
matlab.mixin.Heterogeneous	Superclass for heterogeneous array formation
Enumerations
enumeration	Display class enumeration members and names
meta.EnumeratedValue	Describes enumeration members of MATLAB class
Control Save and Load
save	Save workspace variables to file
load	Load variables from file into workspace
saveobj	Modify save process for object
loadobj	Modify load process for object
Customize MATLAB Behavior
cat	Concatenate arrays along specified dimension
horzcat	Concatenate arrays horizontally
vertcat	Concatenate arrays vertically
empty	Create empty array
disp	Display text or array
display	Display text and numeric expressions
numel	Number of array elements
size	Array dimensions
end	Terminate block of code, or indicate last array index
subsref	Redefine subscripted reference for objects
subsasgn	Subscripted assignment
subsindex	Subscript indexing with object
substruct	Create structure argument for subsasgn or subsref
Custom Object Display
disp	Display text or array
display	Display text and numeric expressions
details	Display array details
matlab.mixin.CustomDisplay	Display customization interface class
matlab.mixin.util.PropertyGroup	Custom property list for object display
Getting Information About Classes and Objects

metaclass	Obtain meta.class object
meta.abstractDetails	Find abstract methods and properties
meta.class.fromName	Return meta.class object associated with named class
meta.package.fromName	Return meta.package object for specified package
meta.package.getAllPackages	Get all top-level packages
meta.class	meta.class class describes MATLAB classes
meta.property	meta.property class describes MATLAB class properties
meta.method	meta.method class describes MATLAB class methods
meta.event	meta.event class describes MATLAB class events
meta.package	meta.package class describes MATLAB packages
meta.DynamicProperty	meta.DynamicProperty class describes dynamic property of MATLAB object
meta.EnumeratedValue	Describes enumeration members of MATLAB class
meta.MetaData	Superclass for MATLAB object metadata
Calling External Functions
Call MEX-File Functions

mexext	Binary MEX-file-name extension
inmem	Names of functions, MEX-files, classes in memory
Call C Shared Libraries

loadlibrary	Load shared library into MATLAB
unloadlibrary	Unload shared library from memory
libisloaded	Determine if shared library is loaded
calllib	Call function in shared library
libfunctions	Return information on functions in shared library
libfunctionsview	Display shared library function signatures in window
libstruct	Convert MATLAB structure to C-style structure for use with shared library
libpointer	Pointer object for use with shared library
lib.pointer	Pointer object compatible with C pointer
Call Java Libraries

javaArray	Construct Java array object
javaclasspath	Return Java class path or specify dynamic path
javaaddpath	Add entries to dynamic Java class path
javarmpath	Remove entries from dynamic Java class path
javachk	Error message based on Java feature support
isjava	Determine if input is Java object
usejava	Determine if Java feature is available
javaMethod	Call Java method
javaMethodEDT	Call Java method from Event Dispatch Thread (EDT)
javaObject	Call Java constructor
javaObjectEDT	Call Java constructor on Event Dispatch Thread (EDT)
cell	Create cell array
class	Determine class of object
clear	Remove items from workspace, freeing up system memory
depfun	List dependencies of function or P-file
exist	Check existence of variable, function, folder, or class
fieldnames	Field names of structure, or public fields of object
im2java	Convert image to Java image
import	Add package or class to current import list
inmem	Names of functions, MEX-files, classes in memory
inspect	Open Property Inspector
isa	Determine if input is object of specified class
methods	Class method names
methodsview	View class methods
which	Locate functions and files
matlab.exception.JavaException	Capture error information for Java exception
Call .NET Libraries

Getting Started
NET.addAssembly	Make .NET assembly visible to MATLAB
NET.isNETSupported	Check for supported Microsoft .NET Framework
NET	Summary of functions in MATLAB .NET interface
enableNETfromNetworkDrive	Enable access to .NET commands from network drive
NET.Assembly	Members of .NET assembly
NET.NetException	Capture error information for .NET exception
Data Types
NET.createArray	Array for nonprimitive .NET types
cell	Create cell array
NET.disableAutoRelease	Lock .NET object representing a RunTime Callable Wrapper (COM Wrapper) so that MATLAB does not release COM object
NET.enableAutoRelease	Unlock .NET object representing a RunTime Callable Wrapper (COM Wrapper) so that MATLAB releases COM object
Properties
NET.setStaticProperty	Static property or field name
Events and Delegates
BeginInvoke	Initiate asynchronous .NET delegate call
EndInvoke	Retrieve result of asynchronous call initiated by .NET System.Delegate BeginInvoke method
Combine	Convenience function for static .NET System.Delegate Combine method
Remove	Convenience function for static .NET System.Delegate Remove method
RemoveAll	Convenience function for static .NET System.Delegate RemoveAll method
Enumerations
bitand	Bit-wise AND
bitor	Bit-wise OR
bitxor	Bit-wise XOR
bitnot	.NET enumeration object bit-wise NOT instance method
Generic Classes
NET.convertArray	Convert numeric MATLAB array to .NET array
NET.createGeneric	Create instance of specialized .NET generic type
NET.invokeGenericMethod	Invoke generic method of object
NET.GenericClass	Represent parameterized generic type definitions
Call COM Objects

actxserver	Create COM server
actxcontrol	Create Microsoft ActiveX control in figure window
actxcontrollist	List currently installed Microsoft ActiveX controls
actxcontrolselect	Create Microsoft ActiveX control from GUI
iscom	Determine whether input is COM or ActiveX object
isprop	Determine whether input is COM object property
get	Get property value from interface, or display properties
set	Set object or interface property to specified value
addproperty	Add custom property to COM object
deleteproperty	Remove custom property from COM object
inspect	Open Property Inspector
propedit	Open built-in property page for control
fieldnames	Field names of structure, or public fields of object
ismethod	Determine whether input is COM object method
methods	Class method names
methodsview	View class methods
invoke	Invoke method on COM object or interface, or display methods
isevent	Determine whether input is COM object event
events	List of events COM object can trigger
eventlisteners	List event handler functions associated with COM object events
registerevent	Associate event handler for COM object event at run time
unregisterallevents	Unregister all event handlers associated with COM object events at run time
unregisterevent	Unregister event handler associated with COM object event at run time
isinterface	Determine whether input is COM interface
interfaces	List custom interfaces exposed by COM server object
release	Release COM interface
delete	Remove COM control or server
move	Move or resize control in parent window
load	Initialize control object from file
save	Serialize control object to file
Call Web Services

callSoapService	Send SOAP message to endpoint
createClassFromWsdl	Create MATLAB class based on WSDL document
createSoapMessage	Create SOAP message to send to server
parseSoapResponse	Convert response string from SOAP server into MATLAB types
Exception Handling
try, catch	Execute statements and catch resulting errors
addCause (MException)	Record additional causes of exception
getReport (MException)	Get error message for exception
last (MException)	Last uncaught exception
rethrow (MException)	Reissue existing exception
throw (MException)	Issue exception and terminate function
throwAsCaller (MException)	Throw exception as if from calling function
MException	Capture error information
Unit Testing Framework
Write Unit Tests

functiontests	Create array of tests from handles to local functions
matlab.unittest.TestCase	Superclass of all matlab.unittest test classes
Run Unit Tests

runtests	Run set of tests
matlab.unittest.TestCase.run	Run TestCase test
matlab.unittest.TestSuite.run	Run TestSuite array using TestRunner object configured for text output
matlab.unittest.TestRunner.run	Run all tests in TestSuite array
matlab.unittest.TestRunner.addPlugin	Add plugin to TestRunner object
matlab.unittest.TestSuite	Class for grouping tests to run
matlab.unittest.Test	Specification of a single test method
matlab.unittest.TestRunner	Class for running tests in matlab.unittest framework
Analyze Test Results

matlab.unittest.TestResult	Result of running test suite
Performance and Memory
Code Performance

bench	MATLAB benchmark
cputime	Elapsed CPU time
memory	Display memory information
profile	Profile execution time for function
profsave	Save profile report in HTML format
tic	Start stopwatch timer
timeit	Measure time required to run function
toc	Read elapsed time from stopwatch
Memory Usage

clear	Remove items from workspace, freeing up system memory
inmem	Names of functions, MEX-files, classes in memory
memory	Display memory information
pack	Consolidate workspace memory
whos	List variables in workspace, with sizes and types
MATLAB Environment Control
commandhistory	Open Command History window, or select it if already open
commandwindow	Open Command Window, or select it if already open
filebrowser	Open Current Folder browser, or select it if already open
workspace	Open Workspace browser to manage workspace
getpref	Preference
setpref	Set preference
addpref	Add preference
rmpref	Remove preference
ispref	Test for existence of preference
matlab.io.saveVariablesToScript	Save workspace variables to MATLAB script
Custom Documentation
builddocsearchdb	Build searchable documentation database
MATLAB API for Other Languages
MATLAB Engine API

mex	Build MEX-function from C/C++ or Fortran source code
MATLAB COM Automation Server

Execute	Execute MATLAB command in Automation server
Feval	Evaluate MATLAB function in Automation server
GetCharArray	Character array from Automation server
PutCharArray	Store character array in Automation server
GetFullMatrix	Matrix from Automation server workspace
PutFullMatrix	Matrix in Automation server workspace
GetVariable	Data from variable in Automation server workspace
GetWorkspaceData	Data from Automation server workspace
PutWorkspaceData	Data in Automation server workspace
MaximizeCommandWindow	Open Automation server window
MinimizeCommandWindow	Minimize size of Automation server window
actxGetRunningServer	Handle to running instance of Automation server
enableservice	Enable, disable, or report status of MATLAB Automation server
Quit	Terminate MATLAB Automation server
MEX-File Creation API

Executable C/C++ MEX-Files
mex	Build MEX-function from C/C++ or Fortran source code
dbmex	Enable MEX-file debugging (on UNIX platforms)
mex.getCompilerConfigurations	Get compiler configuration information for building MEX-files
Call MEX-File Functions
mexext	Binary MEX-file-name extension
inmem	Names of functions, MEX-files, classes in memory
Share MEX-Files
ver	Version information for MathWorks products
computer	Information about computer on which MATLAB software is running
mex.getCompilerConfigurations	Get compiler configuration information for building MEX-files
mexext	Binary MEX-file-name extension
Troubleshoot MEX-Files
dbmex	Enable MEX-file debugging (on UNIX platforms)
inmem	Names of functions, MEX-files, classes in memory
mex	Build MEX-function from C/C++ or Fortran source code
mex.getCompilerConfigurations	Get compiler configuration information for building MEX-files
mexext	Binary MEX-file-name extension
Source Control Integration
checkin	Check files into source control system (UNIX platforms)
checkout	Check files out of source control system (UNIX platforms)
cmopts	Name of source control system
customverctrl	Allow custom source control system (UNIX platforms)
undocheckout	Undo previous checkout from source control system (UNIX platforms)
verctrl	Source control actions (Windows platforms)
Desktop Environment

Startup and Shutdown
matlab (Windows)	Start MATLAB program (Windows platforms)
matlab (UNIX)	Start MATLAB program (UNIX platforms)
exit	Terminate MATLAB program (same as quit)
quit	Terminate MATLAB program
matlabrc	Startup file for MATLAB program
startup	Startup file for user-defined options
finish	Termination file for MATLAB program
Basic Settings
prefdir	Folder containing preferences, history, and layout files
preferences	Open Preferences dialog box
Platform and License
ismac	Determine if version is for Mac OS X platform
ispc	Determine if version is for Windows (PC) platform
isstudent	Determine if version is Student Version
isunix	Determine if version is for UNIX platform
javachk	Error message based on Java feature support
license	Return license number or perform licensing task
usejava	Determine if Java feature is available
ver	Version information for MathWorks products
verLessThan	Compare toolbox version to specified version string
version	Version number for MATLAB and libraries
Help and Support
doc	Reference page in Help browserSearch for term in documentation
help	Help for functions in Command Window
docsearch	Help browser search
lookfor	Search for keyword in all help entries
demo	Access product examples in Help browser
echodemo	Run example script step-by-step in Command Window
