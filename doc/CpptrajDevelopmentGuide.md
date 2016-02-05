CPPTRAJ is code used for processing MD trajectory data as well as other types of data, derived from trajectories or otherwise. CPPTRAJ is a complete rewrite of the PTRAJ code in primarily C++, with the intent being to make the code more readable, leak-free, and thread-safe. The biggest functional change from PTRAJ is the ability to load and process trajectories with different topology files in the same run.

This guide assumes that the reader has at least a basic familiarity with C and C++ object-oriented programming. If you aren’t sure what a constructor is or how pointers work you may have a difficult time coding in Cpptraj. There are several good introduction to C/C++ tutorials on the web that may be helpful.

Coding Conventions
==================

It is important to maintain a consistent coding style within cpptraj so that it remains easy to modify and understand. By following code conventions, it will be easier to read code written by anybody and determine what is happening.

-   Code blocks are indented using 2 spaces. **DO NOT USE TABS** since these are in general not portable between different editors.

-   Try to keep lines between 80 and 100 characters long.

-   Whenever possible, put separate code on separate lines. Exceptions can be made for very simple statements such as logic evaluations and simple initializations. For example,\
    **double x1 = 0.0, x2 = 0.0, x3 = 0.0;**\
    is OK, but\
    **double x1 = var1 <span>\*</span> var2; double x2 = var3 / var4;**\
    is not. There are two reasons: 1) When separate statements share a line it makes using debuggers more difficult, and 2) when separate statements share a line it is harder to read.

-   C++ files have ’.cpp’ suffix, C files have ’.c’ suffix, header files have ’.h’ suffix.

-   All header files should have a ’\#define’ guard to prevent multiple inclusion. The define guard has format:

\#ifndef INC\_\<basefilename\>\_H

\#define INC\_\<basefilename\>\_H

...

\#endif

-   ’using namespace’ should be used sparingly and NEVER in a header file.

-   The order of \#include directives should be (in general): C includes, C++ includes, class definition, any other Cpptraj includes.

-   Use of STL classes/methods is acceptable; use C99 conventions to maximize portability. The only external libraries that should be used are NetCDF and ARPACK/LAPACK/BLAS (both included with AmberTools), i.e. no Boost etc.

-   Do not use iostream for basic IO. All console output should be performed with the functions in CpptrajStdio.h (chiefly mprintf() and mprinterr() for STDOUT and STDERR respectively). All file IO should be performed with CpptrajFile or the derived classes BufferedLine and BufferedFrame. This choice has been made mainly for performance reasons (C file routines are in general much faster than iostream), but also so that all IO is centralized (e.g. CpptrajFile will automatically detect if an input file is compressed). This is also so output can be easily controlled; for example, using mprintf will make sure that during MPI only the master writes.

-   Warnings should be written to STDOUT with mprintf with prefix ’Warning:’; errors should be written to STDERR with mprinterr with prefix ’Error:’.

-   Classes:

    -   Class types are named using *CapWords* (no spaces or underscores, start of each word is a capital letter).

    -   Files containing a class should be named after the class (e.g. ’class TrajectoryFile {};’ in TrajectoryFile.cpp).

    -   Classes which inherit should be named after their base class (e.g. ’class Action\_Distance : public Action { };’).

    -   Public class methods should be listed first; protected methods/variables second; private methods/variables last. All class member variables should be private if possible.

    -   Public class methods are named using *CapWords*.

    -   Private class methods are named using *mixedCase.*

    -   Class variables that are `private` or `protected` are named using *mixedCase\_* (with a trailing underscore).

-   Abbreviations: 1st letter in each word is capitalized. For instance, `Data File List` may be abbreviated `DatFilList` or `DFL.`

-   Variables that have function scope (or lower) and all public variables for classes are named using *mixedCase* (same as *CapWords* except the first letter is lower-case).

-   No one-letter variable names except in loop scopes (e.g. for (int i = 0; i \< N; ++i) { } ), and even then they should be short loops (no more than 10 lines or so).

-   All identifiers in an enumerated type are named using all *CAPS*, and the first identifier should be explicitly initialized (e.g. enum DirectionType { DX = 0, DY, DZ };).

-   There is a *doxygen* rule file to automatically generate code documentation using *doxygen*, so please construct comments in such a doxygen-compatible manner (e.g. JavaDoc etc). See http://www.stack.nl/\~dimitri/doxygen/manual.html for instructions.

Building Cpptraj and Documentation
==================================

Cpptraj is automatically built as part of AmberTools, or it can be built standalone using the configure script in the `$AMBERHOME/AmberTools/src/cpptraj` or `$CPPTRAJHOME` directory. The standalone build is particularly useful for development and testing. Type ’./configure –help’ for a list of configure options. In order to build Cpptraj standalone one needs to specify the location of the NetCDF, zlib, bzlib2, and BLAS/LAPACK/ARPACK libraries if they aren’t in your system path; configure will use the ones in \$AMBERHOME if `’-amberlib’` is specified. The -noX options can be used to disable use of certain libraries.

For example, to build cpptraj standalone:

./configure -amberlib gnu

make install OR cd src && make install

To build the documentation using *doxygen*, you must have *doxygen* installed, and you must have configured AmberTools. Run the command:

make docs

to build the documentation. PDF files and HTML files are generated during this process, showing class inheritance and descriptions from comments written in doxy-format. Open the file `$AMBERHOME/AmberTools/src/cpptraj/doc/html/index.html` to see the class heirarchy and descriptions.

Cpptraj currently lives in 3 key classes:

<span>Cpptraj</span>

:   Defined in main.cpp, controls overall flow (e.g. it is responsible for deciding whether to execute in batch mode or interactive mode).

<span>CpptrajState</span>

:   Defined in **Cpptraj**, holds all of the data, Actions, Analyses, etc.

<span>Command</span>

:   “Static” class used by **Cpptraj** to process user input. The file Command.cpp also contains all of the logic for executing commands.

CpptrajState
============

The main components of **CpptrajState** are:

<span>DataSetList DSL\_;</span>

:   Hold all DataSets. This is essentially how different components can talk to each other, e.g. an Action creates a DataSet in the DataSetList, which can then be used by a subsequent Analysis.

<span>DataFileList DFL\_;</span>

:   Hold all DataFiles. These are either for writing out DataSets or general text output, primarily from Actions/Analyses.

<span>TrajinList trajinList\_;</span>

:   Hold all input trajectories to be processed during a run. Whenever a user inputs a ’trajin’ or ’ensemble’ command, the trajectory/ensemble in question is added to this list. When a ’run’ command is executed, these are the trajectories that are read in a frame at a time so that Actions in the ActionList can process them.

<span>TrajoutList trajoutList\_;</span>

:   Hold output trajectories to be written during a run. This output occurs after all Actions have been processed.

<span>ActionList actionList\_;</span>

:   Hold all Actions to be executed during a run. By default, whenver an Action command is issued the Action in question is initialized and queued up in the ActionList, to be processed during the next run.

<span>AnalysisList analysisList\_;</span>

:   Hold all Analyses to be executed after a run or when a ’runanalysis’ command is given. Similar to ActionList, whenever an Analysis command is issued the Analysis in question is initialized and queued up in the AnalysisList.

Actions
=======

Actions are how Cpptraj derives data from input trajectories. There are two basic classes for Actions:

Action
------

The abstract base class that defines the Action interface. Consists of 4 functions:

<span>Init():</span>

:   Initialize Action, set up DataSets/DataFiles, etc.

<span>Setup():</span>

:   Set up Action for a given Topology.

<span>DoAction():</span>

:   Perform Action on given Frame.

<span>Print():</span>

:   Perform any post-processing or output that occurs outside the main DataSet/DataFile framework.

ActionInit (ActionState.h)
--------------------------

Used to interface with Init(); contains pointers to the master DataSetList and DataFileList in CpptrajState.

<span>DataSetList& DSL(), DataSetList const& DSL() const</span>

:   Reference to the master DataSetList. The appropriate version should be used automatically.

<span>DataSetList<span>\*</span> DslPtr()</span>

:   For Actions that require access to the master DataSetList after Init() (e.g. hbond, which cannot set up hbond time series until hbonds are actually detected in DoAction()), they can store a pointer to the master DataSetList like so:

    (In header): DataSetList<span>\*</span> masterDSL\_;

    (In Action::Init): masterDSL\_ = init.DslPtr();

<span>DataFileList& DFL(), DataFileList const& DFL() const</span>

:   Reference to master DataFileList. The appropriate version should be used automatically.

ActionSetup (ActionState.h) 
----------------------------

Used to interface with Setup(); contains pointers to current Topology and CoordinateInfo, as well as expected number of frames associated with current Topology.

ActionFrame (ActionState.h)
---------------------------

Used to interface with DoAction(); contains pointer to current Frame.

Math-related Classes
====================

Vec3
----

An array of 3 doubles, used to hold XYZ coords. Used for vector math.

Matrix\_3x3
-----------

A 3x3 array of doubles, useful for performing rotations etc. Used for basic matrix math. Can be diagonalized via an internal routine (no need for external math library).

ComplexArray
------------

Used to hold an array of complex numbers. Implemented as a double array instead of using the STL Complex class so that it easily interface with external routines.

PubFFT
------

Interface to FFT routines (either pubfft, which are the FFT routines used by Amber, or FFTW depending on how CPPTRAJ is configured). Currently only 1D forward and backwards FFTs are supported. Makes use of ComplexArray.

Corr.h: CorrF\_Direct, CorrF\_FFT
---------------------------------

Classes used to calculate auto/cross correlation functions from arrays of complex numbers (ComplexArray).

Some Key Classes and Functions
==============================

The following is a brief list of some of the more commonly-used classes and functions in Cpptraj. Classes are more or less self-documented to a certain extent; this section will be focused on how these classes are/should be used.

ArgList
-------

The ArgList class is used throughout Cpptraj. It is the main way that user input is translated to actions, analyses, trajectory IO, etc. Basically, the ArgList class takes a string and separates it into tokens based on a given delimiter or delimiters. For example, the string:

myString = “trajin mytraj.nc 1 100 10”;

can be separated via a space (’ ’) delimeter (the default) into 5 tokens like so:

ArgList myArgs(myString);

The resulting ArgList internally looks something like:

0: trajin

1: mytraj.nc

2: 1

3: 100

4: 10

A custom delimeter string containing 1 or more characters can also be used. For example, the following string:

myString = “d01,d02,d03,d04”;

can be separated via a comma (’,’) delimiter into 4 tokens like so:

ArgList myArgs(myString, “,”);

The resulting ArgList internall lookws something like:

0: d01

1: d02

2: d03

3: d04

These tokens (or arguments) are stored internally as an STL vector of strings. ArgList provides many functions to access user arguments. A second array of boolean values records whether an argument has been accessed. This concept is functionally similar to the argumentStack in Ptraj; however, it avoids the constant memory allocation/deallocation when arguments are added/accessed, and allows an argument list to be re-used if desired. The two main ways arguments are usually accessed are through “GetNextX” and “GetKeyX” functions.

“GetNext” functions return the next argument of the desired type. For example, using the ArgList created in the first example from “trajin mytraj.nc 1 100 10” and assuming all arguments are unmarked, GetStringNext() would return “trajin”, while getNextInteger() would return “1”; in both cases the argument returned would be marked, so that a subsequent call to GetStringNext() would return “mytraj.nc” and so on. Another very commonly used “GetNext” function is the “GetMaskNext()” function, which returns the next atom mask expression (so noted because it will begin with ’:’, ’@’, ’<span>\*</span>’); an example of this will be shown below.

“GetKey” functions return an argument next to a specified “key” string. Take for example the argument list created from “rmsd R1 @CA ref <span>[</span>myref<span>]</span> out rmsd.dat”:

0: rmsd

1: R1

2: @CA

2: ref

3: <span>[</span>myref<span>]</span>

4: out

5: rmsd.dat

If we want to access a specific argument, we use a “GetKey” function. For example, if we want to know the filename specified by ’out’, we would use GetStringKey(“out”); this would return “rmsd.dat”, and mark both “out” and “rmsd.dat”. Similarly, GetStringKey(“ref”) would return “<span>[</span>myref<span>]</span>”. At this point we could also use the GetMaskNext() function to get the atom mask expression “@CA”.

### ArgList Example

ArgList myArgs(“test\_command cutoff 2.0 nval 3 name MyTest :2-30@CA extra”);

double Cut = myArgs.getKeyDouble(“cutoff”, 0.0); // Value 2.0

int Nval = myArgs.getKeyInteger(“nval”, 0);      // Value 3

int Ntypes = myArgs.getKeyInteger(“ntypes”, 0);  // Value 0

std::string Name = myArgs.GetStringKey(“name”);  // Value “MyTest”

std::string Out = myArgs.GetStringKey(“out”);    // Empty

std::string maskExp = myArgs.GetMaskNext();      // Value “:2-30@CA”

std::string mask2Exp = myArgs.GetMaskNext();     // Empty

// At this point only “extra” will be unmarked.

Topology
--------

The Topology class describes how a system is laid out in terms of Atoms, Residues, and Molecules (all of which are classes themselves). It may also hold parameters which describe interactions between Atoms (e.g. bonds, angles, dihedrals, etc). The Topology class is chiefly used in Trajectory input/output and setting up atom masks (see below).

The Topology class has several routines that return strings of atom and residue names:

<span>std::string TruncResAtomName(int atom)</span>

:   Format: “\<res name\>\<res num\>@\<atom name\>”

<span>std::string AtomMaskName(int atom)</span>

:   Format: “:\<res num\>@\<atom name\>”

<span>std::string TruncAtomNameNum(int atom)</span>

:   Format: “\<atom name\>\_\<atom num\>”

<span>std::string TruncResNameNum(int residue)</span>

:   Format: “\<res name\>:\<res num\>”

### Examples

Iterate over all atoms in a certain residue. This can be accomplished like so:

// Iterate over all atoms in residue 4, print charge.

for (int atom\_index = Top.Res(4).FirstAtom;

         atom\_index != Top.Res(4).LastAtom;

       ++atom\_index)

  mprintf(“Atom %i charge= %g\\n”, atom\_index+1, Top<span>[</span>atom\_index<span>]</span>.Charge());

Iterate over all atoms bonded to a certain atom and pick out the hydrogens:

// Iterate over all atoms bonded to atom 66

for (Atom::bond\_iterator bond\_atom = Top<span>[</span>66<span>]</span>.bondbegin();

                         bond\_atom != Top<span>[</span>66<span>]</span>.bondend();

                       ++bond\_atom)

  if (Top<span>[</span><span>\*</span>bond\_atom<span>]</span>.Element() == Atom::HYDROGEN)

    mprintf(“Hydrogen %s bonded to atom %s\\n”,

            Top.AtomMaskName(<span>\*</span>bond\_atom).c\_str(),

            Top.AtomMaskName(66).c\_str());

AtomMask/CharMask
-----------------

The AtomMask and CharMask classes keep track of what atoms for a given Topology are selected based on a given mask expression. The AtomMask class holds information on selected atoms only, while CharMask has the state of all atoms (selected/not selected). AtomMask is as an integer mask, where the atom numbers currently selected are stored as an array of integers. Since one is usually only interested in selected atoms, most times AtomMask is all that is needed and so is the most used mask class in Cpptraj; for example, all routines that take masks in the Frame class use the AtomMask class.

The CharMask class has an internal character array where the state of each atom is stored. It has functions that can then be used to interrogate if a certain atom or atoms are selected or not.

In typical use, there are 3 phases to using AtomMask or CharMask: 1) initialization with a mask expression, 2) setup via a Topology class, and 3) iteration over the mask/interrogation of the mask. Initialization with a mask expression performs all necessary tokenization of the mask expression string and prepares the mask to be set up, but does not actually select atoms. The mask expression can be used during AtomMask construction or passed in via SetMaskString():

AtomMask<span>\*</span> Mask = new AtomMask(“@CA”);

AtomMask Mask(“@CA”);

AtomMask Mask; Mask.SetMaskString(“@CA”);

Setup occurs via a Topology class (since in order to set up a mask you need to know atom names/numbers, residue name/number/types etc). This can be done using SetupIntegerMask() or SetupCharMask() to set up an integer mask (more common) or a char mask:

AtomMask iMask(“@CA”);

Top.SetupIntegerMask( iMask );

CharMask cMask(“:1-20”);

Top.SetupCharMask( Mask );

Once a mask has been setup the Nselected() function returns the number of selected atoms, while the None() function returns true if no atoms were selected. If necessary, one can convert between the mask types post-setup by using AtomMask::ConvertToCharMask() / CharMask::ConvertToIntMask() routines:

AtomMask mask( CharMask.ConvertToIntMask(), CharMask.Natom() )

CharMask mask( AtomMask.ConvertToCharMask(), AtomMask.Nselected() )

The final stage is to make use of the atom mask. One can iterate over selected atoms in an integer atom mask using the STL-like const\_iterator variable and begin() and end() functions - this is the recommended way to use atom masks:

for (AtomMask::const\_iterator atomnum = Mask.begin();

                              atomnum != Mask.end(); 

                            ++atomnum)

  mprintf(“Selected atom %i\\n”, <span>\*</span>atomnum);

One can also access members of the integer array directly via the bracket (’<span>[</span><span>]</span>’) operator:

for (int maskidx = 0; maskidx \< Mask.Nselected(); ++maskidx)

  mprintf(“Selected atom %i\\n”, Mask<span>[</span>atomidx<span>]</span>);

To see if atom(s) are selected in a CharMask, use the AtomInCharMask() and AtomsInCharMask() functions. The former returns true if a specified atom is selected, the latter returns true if any atoms within a given range are selected. For example:

for (int atom = 0; atom \< Top.Natom(); ++atom)

  if (Mask.AtomInCharMask(atom)) 

    mprintf(“Selected Atom %i\\n”, atom);

for (int rnum = 0; rnum \< Top.Nres(); ++res)

  if (Mask.AtomsInCharMask( Top.Res(rnum).FirstAtom(),

                            Top.Res(rnum).LastAtom() ))

    mprintf(“Selected Residue %i\\n”, rnum);

Frame
-----

The Frame class is in many ways the workhorse of Cpptraj, as it holds all XYZ coordinates for a given input frame, and optionally box coordinates, masses, replica indices, temperature, time, and/or velocities. Note that although mass is stored in Topology, it is also stored in Frame since many calculations require it (center of mass, mass-weighted RMSD, etc). Coordinates and velocities are stored with double precision. Many routines are available to do things like calculate the center of mass of atoms, rotate, translate, scale, and so on. A major use of the Frame class is to perform RMSD calculations.

Frames are typically set up in two phases. The first phase is memory allocation, which occurs via constructors or the SetupFrameX routines. This should be done as little as possible since memory allocation is relatively expensive. The second phase is actually setting the coordinates, which occurs via the SetX routines. For example, the following code will set up a Frame (newFrame) with the coordinates from another Frame (oldFrame) based on a previously set up AtomMask (mask) and Topology (top) corresponding to oldFrame:

Frame newFrame;

// Allocate memory, copy in masses based on mask.

newFrame.SetupFrameFromMask( mask, top.Atoms() );

// Set only coordinates from oldFrame based on mask.

newFrame.SetCoordinates( oldFrame, mask );

Alternatively, this can be done in one step with a constructor:

Frame newFrame( oldFrame, mask );

The advantage of separating out Setup and Set is that memory reallocation is kept to a minimum. For example, if we wanted to use newFrame to hold a different set of coordinates (of the same size as newFrame or smaller) we might do something like:

newFrame.SetCoordinates( differentFrame );

If the new coordinates might be bigger than the current size of newFrame, we could explicitly call a SetupFrameX routine; this will only reallocate if the new size is greater than the current maximum size. A Frame remembers the largest size it was ever allocated for, so reallocation is kept to a minimum.

There are two basic ways to access coordinates within Frame:

<span>const double<span>\*</span> XYZ(atom)</span>

:   return pointer to beginning of XYZ coordinates for given atom (max Natom()).

<span>const double<span>\*</span> CRD(coord)</span>

:   return pointer to given coordinate (max size()).

Note that you can get pointers to the raw coordinates, but it is not recommened to use these in general.

### Using Frame for RMSD calculations

Unlike some of the other functions of Frame, the RMSD functions do not take a mask - it is assumed all atoms in the Frame are involved in the RMSD calculation. This is done for performance reasons. If a subset of atoms is desired for an RMSD calculation the reference and target Frames should be modified beforehand. Since the reference structure usually does not change it is often beneficial to pre-center the reference at the origin. For example, given a reference Frame (Ref), a Farget frame (Tgt), and an AtomMask (mask):

bool useMass = false;

top.SetupIntegerMask( mask );

Frame selectedRef, selectedTgt;

// Set up and pre-center reference.

selectedRef.SetupFrameFromMask( mask, top.Atoms() );

selectedRef.SetCoordinates( Ref, mask );

// refTrans will contain translation from origin to reference.

Vec3 refTrans = selectedRef.CenterOnOrigin( useMass );

// Set up target.

selectedTgt.SetupFrameFromMask( mask, top.Atoms() );

selectedTgt.SetCoordinates( Tgt, mask );

// Calculate RMSD. tgtTrans is translation from target to origin.

Vec3 tgtTrans;

Matrix\_3x3 rot\_matrix;

double rmsd = selectedTgt.RMSD\_CenteredRef( selectedRef,

                                            rot\_matrix,

                                            tgtTrans,

                                            useMass );

// Best-fit rotate/translate current Target to Reference.

Tgt.Trans\_Rot\_Trans( tgtTrans, rot\_matrix, refTrans );

Now for subsequent RMS calculations to the same Reference, only the selected Target frame needs to have its coordinates set:

selectedTgt.SetCoordinates( Tgt2, mask );

rmsd = selectedTgt.RMSD\_CenteredRef( selectedRef, ...

Console and File Input/Output
=============================

Output to STDOUT/STDERR: CpptrajStdio.h
---------------------------------------

The file CpptrajStdio.cpp contains functions used to write output to standard output (STDOUT) and standard error (STDERR). This is accomplished with the C printf-like functions mprintf() and mprinterr() respectively:

mprintf(const char<span>\*</span> format, ...); // STDOUT

mprinterr(const char<span>\*</span> format, ...); // STDERR

The ’m’ prefix stands for “master”, and ensures that when CPPTRAJ is running via MPI that only the master thread is able to write with these functions (note the same does NOT apply for OpenMP). The syntax is the same as basic printf - a format string followed by any variables. There are numerous resources that describe printf syntax in detail. Some useful syntax is listed here:

<span>%i</span>

:   Integer

<span>%f</span>

:   Floating point number

<span>%g</span>

:   Use scientific or floating point representation, whichever is shorter.

<span>%s</span>

:   Character string

<span>%c</span>

:   Single character

<span>\\n</span>

:   Newline

<span>\\t</span>

:   Tab

For example:

\#include “CpptrajStdio.h”

int myInteger = 3;

double myDouble = 5.43;

string myString = “easy”;

mprintf(“Using printf is %s; %i is an integer and %f is a double.\\n”,

        myString.c\_str(), myInteger, myDouble);

Note that the string function c\_str() must be used to print C++ strings. Output:

Using printf is easy; 3 is an integer and 5.430000 is a double.

CpptrajFile
-----------

The CpptrajFile class provides basic file input and output operations. It can handle reading and writing both Gzip and Bzip2 compressed files, and through the FileName class performs tilde-expansion on file names (via globbing) as well as separates the file name into its base name, extension, and compressed extension. The file can be opened immediately, or set up first and then opened later. Once it has been set up it can be opened or closed multiple times. The CpptrajFile class destructor will automatically close the file if it is open at time of destruction.

BufferedLine
------------

The BufferedLine class is a child of CpptrajFile used for text files that will be read in line by line (note that writing is not possible with this class). The class has an internal buffer, which chunks of the input file are read into. The Line() routine can be used to read that chunk line by line; this avoids potentially expensive file IO. When the chunk is empty a new chunk is read in. The line can be further split into Tokens (similar to ArgList) and read one token at a time; this can be useful for e.g. determining the number of columns in a file.

BufferedFrame
-------------

The BufferedFrame class is a child of CpptrajFile used for highly-formatted text files that will be read/written multiple lines at a time (such as the Amber ASCII trajectory format). This class is set up for a certain total number of elements of a certain character width with a certain number of elements per line, which can then be read to or written from a character buffer in one entire chunk.

FileName, FileName.h
--------------------

The FileName class is used by all file-related classes. It is more powerful than a string and will automatically do tilde expansion and split the name into path, base name, extension, etc. FileName.h also contains the File namespace which includes File::Exists() for testing whether a file can be opened and File::NameArray and File::ExpandToFilenames() for getting an array of files using wildcard matching.

Trajectory Input/Output
=======================

Trajectory input and output (IO) is handled via high-level and low-level classes. At the highest level trajectory input is provided by Trajin-derived classes (Trajin\_Single and Trajin\_Multi) for reading one frame at a time (i.e. during ’trajin’ runs) and EnsembleIn-derived classes (EnsembleIn\_Single, which is experimental currently, and EnsembleIn\_Multi) for reading multiple frames at a time (i.e. during ’ensemble’ runs). Trajectory output is currently handled by Trajout\_Single for writing one frame at a time to a single file, and EnsembleOut-dervived classes for writing multiple frames out at a time. At the lower level IO is handled by format-specific classes which inherit from TrajectoryIO (TrajectoryIO.h), which are called Traj\_X by convention (e.g. Traj\_AmberNetcdf for Amber NetCDF trajectories). TrajectoryIO classes are contained and set up within the higher level classes. Note that there is currently no Trajout\_Multi (write a single frame to multiple files) since this functionality is already handled by TrajoutList.

Trajin\_Single
--------------

This is for reading in a single frame at a time from a single file.

Trajin\_Multi
-------------

This is for reading in a single frame at a time from multiple files (e.g. getting a frame at a specifed temperature from a T-REMD ensemble).

EnsembleIn\_Single
------------------

For reading in multiple frames at a time from a single file. Currently experimental.

EnsembleIn\_Multi
-----------------

For reading in multiple frames at a time from multiple files, optionally sort the frames.

Trajout\_Single
---------------

For writing out a single frame at a time to a single file.

EnsembleOut\_Single
-------------------

For writing out multiple frames at a time to a single file. Currently experimental.

EnsembleOut\_Multi
------------------

For writing out multiple frames at a time to multiple files. Used for ’ensemble’ run trajectory output and LES trajectory splitting.

Topology Input/Output
=====================

The *ParmIO* class is a base class for all topology file formats. This provides an easy mechanism for extracting the system topology from any number of file formats. The *ParmFile* class is a wrapper around the *ParmIO* classes that hides the implementation details for each data file type from you. You should interact with *ParmIO* objects through *ParmFile* handlers. *ParmFile* provides the ability to both read and write topology file objects of any class.

One thing that sets *ParmIO* and *ParmFile* apart from *TrajectoryIO/TrajectoryFile* and *DataIO/DataFile* is its connection with the *Topology* class. *Topology* objects contain as much of the information in the Amber topology file as can be parsed from the information present in the ParmIO object (and figured out based on atomic arrangements). A *Topology* instance is the first argument passed to the *ParmFile::Read* function, followed by the name of the topology file. Unlike the *DataFile* and *TrajectoryFile* classes, *ParmFile* does not have a reference to the *ParmIO* object to forward read/write information to. It exists simply to fill the *Topology* class with the relevant data structures and inform it how to do the rest. The *Topology* class is format-independent, providing a layer of abstraction to make other parts of the code that require topology information less error-prone while coding.

Every *ParmIO* subclass implements a *ReadParm* method that takes a *Topology* instance as the first argument and fills as much of the information there as possible. Afterwards, the CommonSetup method of the Topology class is called to finish setup and determine bond information (from atom distances if not present directly in the file format) and molecule information (based on the bonded structure). The currently available types of topologies are summarized in .

<span>|c|c|</span> *ParmIO* Subclass (*ParmFormatType*) & Description<span>\
</span>Parm\_Amber (AMBERPARM) & Amber style topology file (OLD and NEW styles)<span>\
</span>Parm\_CharmmPsf (CHARMMPSF) & CHARMM PSF topology file format (used by NAMD, too)<span>\
</span>Parm\_Mol2 (MOL2FILE) & TRIPOS Mol2 file<span>\
</span>Parm\_PDB & PDB File<span>\
</span>

[tbl:Cpptraj-Parm-Formats]

The DataSet and DataFile Framework
==================================

One of the goals in writing cpptraj was to try and generalize data collection and output, so that any action could output any generated data in any format known to cpptraj without having to write any extra code. For example, data generated by a distance calculation can be output in columns, as a Grace file, as a Gnuplot file, or as all three. To that end Actions and Analyses have access to the main DataSetList in CpptrajState (usually named DSL inside Actions) and the main DataFileList (usually named DFL inside actions). A DataSet can be generated by an Action, but because it is held outside the Action in the master DataSetList it can persist after the Action that generated it has been destroyed and be used in subsequent Analyses etc. DataSets have a base type which determines what kind of data it can hold, and can also belong to a group; DataSets in the same group behave in a similar fashion. For example, all DataSets in the COORDINATES group hold coordinates, etc.

The MetaData class
------------------

DataSets are associated with meta data which is used to both describe and categorize the data. DataSet selection uses meta data; meta data is also used to sort DataSets. This information is contained within the MetaData class. A DataSet can have a name, aspect, index, type etc. This allows for DataSets to e.g. be marked as an alpha torsion (M\_TORSION, ALPHA). There are currently 9 MetaData variables:

<span>name\_</span>

:   The DataSet name; this is the most general level of classification. Used in searches.

<span>fileName\_</span>

:   If the data was read in from a file, the name of that file. Used in searches; both the full path and base name can be used to match (e.g. for reference frame data sets).

<span>aspect\_</span>

:   A “sub-name”, used to differentiate between different aspects of similar data. For example, DataSets generated by the ’nastruct’ Action have different aspects for nucleic acid “stretch”, “shear”, “stagger”, etc. Used in searches.

<span>legend\_</span>

:   This is the name that will be used when writing data out to a file. Not used in searches. A default one is created if one is not provided.

<span>idx\_</span>

:   (Index) A number which can be used to further differentiate data. It could correspond to residue, atom, etc. Used in searches.

<span>ensembleNum\_</span>

:   (Ensemble number) For use during ensemble processing; this is set by the DataSetList to differentiate data from different ensembles.

<span>scalarmode\_</span>

:   Internal categorization of the type of data, e.g. a distance, angle, torsion, etc. Certain functions will use this data - for example data sets marked as angle, torsion, or pucker will take periodicity into account when averaging (DataSet\_1D::Avg()).

<span>scalartype\_</span>

:   Internal sub-type; so in the case of a torison what type of torsion it is (phi, psi etc) or in the case of a matrix what kind of matrix, etc.

<span>timeSeries\_</span>

:   Whether the DataSet is a time series; used by DataSetList to determine whether to call the Allocate() function.

MetaData can be set using the DataSet::SetMeta() routine, but it is recommended that this be done sparingly once a DataSet is part of the master DataSetList since this could create conflicts.

The TextFormat class
--------------------

This class is used when writing the data set out to text files. It creates a printf-like format string of a given type with specified width and precision.

Brief DataSet/DataFile Example
------------------------------

A simple usage example is given here. Say for example we want to track a distance in an Action. The first step is to create the DataSet in the Init() routine.

dist\_ = DSL-\>AddSet(DataSet::DOUBLE, MetaData(actionArgs.GetStringNext(),

                                              MetaData::M\_DISTANCE), Dis);

if (dist\_==0) return 1;

In the first line a DataSet class of type DOUBLE is added to the master DataSetList. The various types are enumerated in DataSet::DataType (DataSet.h). The DataSet will be named whatever the next string is in the actionArgs ArgList. If there is no name, a default one will be created based on the given default Dis and the DataSet’s overall position in the DataSetList (so in this case the default could be something like Dis\_00000). The DataSet is also given the scalarMode M\_DISTANCE, which is information that other Actions/Analyses can use (like Analysis\_Statistics). What is returned is a pointer to the DataSet; DataSet is actually a base class that specific DataSet types inherit (in this case DataSet\_double). In this way the interface is generalized.

The next step is to add the DataSet to the DataFileList.

DataFile<span>\*</span> outfile = DFL-\>AddFile(distanceFile, actionArgs);

if (outfile != 0) outfile-\>AddDataSet( dist\_ );

In the first line a pointer to a new or exsiting (if already created somewhere else) DataFile is returned only if the string distanceFile is not empty; this allows specification of output files to be optional. The actionsArgs ArgList is passed in so that the DataFile outfile can process any DataFile-related arguments. In the second line the DataSet is added to the DataFile. In this way output from multiple actions can be combined rather than overwritten. The machinery of the DataFileList takes care of output (formatting etc) from there.

The final phase is actually adding data to the DataSet. So for example in the action() routine you could have:

double distance = sqrt( DIST2\_NoImage( V1, V2 ) );

dist\_-\>Add(frameNum, &distance);

In the first line the value ’distance’ is being calculated. In the next line the value from ’distance’ is being added to DataSet dist\_ with frame number ’frameNum’ (automatically set within Action). Notice that the address of ’distance’ is passed rather than the value; this is a necessity from the generalization of the DataSet interface. DataSet has no idea a prior what the data type might be, so in the Add routine the value is cast to what the underlying DataSet implementation expects. This allows the Add routine to be used for double, float, int, string, etc. This DataSet could also be cast back to it’s actual type to access other routines, e.g.:

DataSet\_double& ds\_dist = static\_cast\<DataSet\_double&\>( <span>\*</span>dist\_ );

DataSet\_1D / SCALAR\_1D
------------------------

Base class for 1D scalar data sets (like DataSet\_double, DataSet\_Mesh, etc). Basically hold a series of numbers.

DataSet\_2D / MATRIX\_2D
------------------------

Base class for 2D matrices.

DataSet\_3D / GRID\_3D
----------------------

Base class for 3D grids.

DataSet\_Coords / COORDINATES
-----------------------------

Base class for DataSets which hold coordinates.

Most development for Cpptraj will likely be in adding new functionality; actions, analyses, and trajectory/topology/data file formats. This part of the manual will provide guidance and some helpful hints to this end. In general, adding new functionality is done by writing an implementation of the desired class type (e.g. for actions, inherit from the Action class) and then adding that class to the container for that specific functionality (e,g, in the case of actions, ActionList).

Adding Actions - Example
========================

All actions inherit from the Action abstract base class. The Action class itself inherits from the DispatchObject class so that it can be associated with an allocator (to create the action) and a help function. There are four functions that every action must implement: Init(), Setup(), DoAction(), and Print(). Init() is called when the action is first created, and processes input arguments, sets up DataSets/DataFiles, deals with reference frames, and sets the debug level. Setup() is called to set the action up for a specific topology, and so handles anything Topology-related (such as parsing atom masks). The DoAction() function is called to actually perform the action on input coordinate frames. The Init(), Setup(), and DoAction() functions return a special type of integer, Action::RetType, which described the result of the action:

<span>Action::OK</span>

:   Action is successful.

<span>Action::ERR</span>

:   Action is not successful.

<span>Action::USE\_ORIGINAL\_FRAME</span>

:   Action requests that the original unmodified topology/frame be used (see e.g. Action\_Unstrip in Action\_Strip.h).

<span>Action::SUPPRESS\_COORD\_OUTPUT</span>

:   Action requests that further processing of the current coordinate frame be skipped (see e.g. Action\_RunningAvg).

<span>Action::SKIP</span>

:   Non-fatal problem occurred during setup; skip Action until next Setup call.

<span>Action::MODIFY\_TOPOLOGY</span>

:   Setup routine has modified the Topology/CoordinateInfo.

<span>Action::MODIFY\_COORDS</span>

:   DoAction routine has modified the Frame.

The final function is Print(), which is called after all trajectory processing is complete and performs any additional calculation or output necessary. This function can be blank if such functionality is not needed, but it still must be implemented.

In addition to Action, there are currently two additional action-related classes that actions may want to inherit from. The ImagedAction class is for classes that may need to calculate imaged distances, and the ActionFrameCounter class is for actions that may want to process subsets of input frames (see e.g. the Action\_Matrix action).

As an example, we will go through the creation of a simplified version of the Action\_Distance class for calculating distances; this will cover using the DataSet, DataFile, AtomMask, and ImagedAction classes as well.

Create the Class Header
-----------------------

As mentioned in the style guide, header files should be named after the class, so the Action\_Distance class will go in a file named “Action\_Distance.h”. The first thing to do is create a “header guard” - this will prevent issues with multiple inclusion. The header guard should be named after the class and header file, so for Action\_Distance.h:

\#ifndef INC\_ACTION\_DISTANCE\_H

\#define INC\_ACTION\_DISTANCE\_H

Next comes the class description. Since distance calculations may involve imaging we also include the ImagedAction class as a variable to simplify image handling:

class Action\_Distance: public Action {

Following the style guide, we first implement any public methods. For actions this is at least the constructor, the allocator (named Alloc() by convention), and the Help() function. The allocator and help functions need to be static so that they can be called without instantiating the class.

public:

  Action\_Distance();  ///\< Constructor

  /// Allocator

  static DispatchObject<span>\*</span> Alloc() { return (DispatchObject<span>\*</span>)new Action\_Distance(); }

  static void Help(); ///\< Help function

The implemented functions Init(), Setup(), DoAction(), and Print() can be either public or private, although the preference is private.

The private section is where all functions and variables specific to the class will go. First, we add entries for the functions inherited from the Action base class which must be implemented. Since we will not need to do any post-processing for this action, the Print() function is empty:

  Action::RetType Init(ArgList&, ActionInit&, int);

  Action::RetType Setup(ActionSetup&);

  Action::RetType DoAction(int, ActionFrame&);

  void Print() {}

Next we define the class variables. For Action\_Distance we will want a DataSet to hold the calculated distances, two AtomMasks to describe the points between which the distance should be calculated, and a variable to indicate whether the distance should be mass-weighted.

private:

  DataSet<span>\*</span> dist\_;      ///\< Will hold DataSet of calculated distances.

  bool useMass\_;       ///\< If true, mass-weight distances.

  AtomMask Mask1\_;     ///\< First atom selection.

  AtomMask Mask2\_;     ///\< Second atom selection

  ImagedAction image\_; ///\< Holds imaging info

All variables related to imaging are already include via the ImagedAction class.

Last, end the class definition and finish the header guard:

};

\#endif

Create the Class Implementation
-------------------------------

Following the naming scheme, the class implementation will go into Action\_Distance.cpp. The first part of this file will have the necessary \#include directives. We need \<cmath\> for the square root function, Action\_Distance.h for the class definition, and CpptrajStdio.h for wiriting to the console.

\#include \<cmath\>

\#include Action\_Distance.h

\#include CpptrajStdio.h

First we will need to create the class constructor. It is encouraged that users make use of initalizer lists (which tend to be more efficient) for this purpose. In this case we have two non-class variables: dist\_, which is a pointer to a DataSet, and useMass\_, which is boolean:

Action\_Distance::Action\_Distance() : dist\_(0), useMass\_(true) {}

Next, ensure that the Help() function has an implementation. Note that in cpptraj “mprintf” is used over “printf” for making any future IO modifications easier:

void Action\_Distance::Help() {

  mprintf(distance <span>[</span>\<name\><span>]</span> \<mask1\> \<mask2\> <span>[</span>out \<filename\><span>]</span> <span>[</span>geom<span>]</span> <span>[</span>noimage<span>]</span>\\n);

}

### Init() - Parse user arguments, set up DataSets/DataFiles etc

Init() is called when the action is created and is responsible for parsing the Argument list (ArgList) and inital setup. Init() has the same input arguments for every action:

Action::RetType 

  Action\_Distance::Init(ArgList& actionArgs, ActionInit& init, int debugIn)

{

The input arguments are as follows: **actionArgs** (ArgList class) contains arguments from user input, **init** (ActionInit class) contains the master DataSetList and master DataFileList, and **debugIn** is the current debug level for all actions. It is up to the action implementation whether it wants to record the debug level or not.

Typical order of argument processing is keywords, masks, DataSet name. First we will process the keywords ’noimage’, ’geom’, and ’out \<filename\>’.

In order to determine whether the action will try to use imaging we call the InitImaging() function (inherited from the ImagedAction class). If the ArgList actionArgs contains the string “noimage”, false will be sent to InitImaging to disable imaging:

  image\_.InitImaging( !(actionArgs.hasKey(noimage)) );

Next we set useMass\_. If actionArgs contains the string “geom”, useMass\_ will be set to false:

  useMass\_ = !(actionArgs.hasKey(geom));

Next, we will try to create an output DataFile:

  DataFile<span>\*</span> outfile = init.DFL().AddDataFile( actionArgs.GetStringKey(out),

                                              actionArgs );

The behavior of AddDataFile() depends on the result from actionArgs.GetStringKey(); if “out” is present in actionArgs, the next string (presumably \<filename\>) is returned and passed to AddDataFile(), and a DataFile will be returned corresponding to \<filename\>. If “out” is not present nothing will be returned, no file will be set up, and outfile will be null (0).

Next, we will get two atom mask expressions. We will require that the user must specify two masks, so if either of the strings is empty return an error:

  std::string mask1 = actionArgs.GetMaskNext();

  std::string mask2 = actionArgs.GetMaskNext();

  if (mask1.empty() || mask2.empty()) {

    mprinterr(Error: distance: Requires 2 masks\\n);

    return Action::ERR;

  }

Now we can use the mask expression strings to initialize the two AtomMask classes (note that this tokenizes the mask expressions but does not yet set them up since we need topology information to do that):

  Mask1\_.SetMaskString(mask1);

  Mask2\_.SetMaskString(mask2);

Next we will use the master DataSetList (init.DSL()) to create a DataSet to store the calculated distances. We will use a version of DataSetList::AddSet() that allows us to specify the DataSet type, MetaData, and a default name if no name is specified (“Dis”). If any errors occur in creating the DataSet, NULL (0) will be returned. Note that the string returned by actionArgs.GetStringNext() is implicitly converted to a MetaData class.

  dist\_ = init.DSL().AddSet(DataSet::DOUBLE, actionArgs.GetStringNext(), Dis);

  if (dist\_==0) return Action::ERR;

If a DataFile was previously set up, we now add the DataSet to this DataFile:

  if (outfile != 0) outfile-\>AddDataSet( dist\_ );

Last, we print out some information regarding how the Action has been initialized and return Action::OK to indicate successful intialization:

  mprintf(    DISTANCE: %s to %s,Mask1\_.MaskString(), Mask2\_.MaskString());

  if (!image\_.UseImage())

    mprintf(, non-imaged);

  if (useMass\_)

    mprintf(, center of mass);

  else 

   mprintf(, geometric center);   mprintf(.\\n);

  return Action::OK;

}

**IMPORTANT:** Note that this is the only time in which the master DataSetList is passed to the Action. If the Action will need to set up DataSets later (because e.g. they may depend on what’s in the Topology, like in the case of the ***multidihedral*** command), it should save a pointer to the master DataSetList using the **init.DslPtr()** function, e.g.

masterDSL\_ = init.DslPtr();

### Setup() - Set up Topology-related parts of the Action

Setup() is called whenever the action needs to be set up for a given Topology file. Any component of the action that depends on Topology (in this case the AtomMasks and the Imaging) is handled here. The arguments to Setup() are:

Action::RetType Action\_Distance::Setup(ActionSetup& setup) {

Note that the **setup** variable (ActionSetup class) contains a pointer to the current Topology and current trajectory CoordinateInfo, as well as the number of expected frames associated with this Topology during the current run. Actions that want to modify the current Topology or CoordinateInfo for subsequent Actions can do so using the **setup** variable (see e.g. Action\_Strip).

First, we setup the AtomMasks. Each AtomMask is passed to the current topology using the SetupIntegerMask() function, which will create an integer array containing only the selected atoms based on the mask expression. If we needed to know both selected and unselected atoms we could use the SetupCharMask() function instead.

  if (setup.Top().SetupIntegerMask( Mask1\_ )) return Action::ERR;

  if (setup.Top().SetupIntegerMask( Mask2\_ )) return Action::ERR;

After this, we print some information about what atoms are selected (note we could also use the MaskInfo() function of AtomMask for this). For calculating distance, we need to make sure atoms were actually selected (using the None() function of AtomMask). If no atoms were selected this may be because the mask is only valid for certain Topologies during the run, so in that case make it a non-fatal error (i.e. a Warning) and return Action::SKIP:

  mprintf(\\t%s (%i atoms) to %s (%i atoms),Mask1\_.MaskString(), Mask1\_.Nselected(),

          Mask2\_.MaskString(),Mask2\_.Nselected());

  if (Mask1\_.None() || Mask2\_.None()) {

    mprintf(\\nWarning: distance: One or both masks have no atoms.\\n);

    return Action::SKIP;

  }

Next we determine if imaging can actually be performed based on the box information present in the current trajectory’s CoordinateInfo; if there is no box information imaging cannot be performed. We do this with the image\_.SetupImaging() function (ImagedAction class). The image\_.ImagingEnabled() function will let us know if imaging for this Topology is possible or not:

  image\_.SetupImaging( setup.CoordInfo().TrajBox().Type() );

  if (image\_.ImagingEnabled())

    mprintf(, imaged);

  else

    mprintf(, imaging off);

  mprintf(.\\n);

Now all Topology-dependent aspects of the action are set up. Return Action::OK.

  return Action::OK;

}

**IMPORTANT:** Note that this is the only time in which a Topology is passed to the Action. If the Action requires Topology information later (such as in DoAction() or Print()) it should save a pointer to the Topology using the **setup.TopAddress()** function, e.g.

currentParm\_ = setup.TopAddress();

### DoAction() - Process input Frame

Coordinates are read in a frame at a time and stored in a Frame class, which is then passed to each action in the ActionList. The DoAction() function is called to process a coordinate Frame. The arguments are:

Action::RetType Action\_Distance::DoAction(int frameNum, ActionFrame& frm) {

The first argument **frameNum** is the current frame number (starting at 0). Note that the **frm** variable (ActionFrame class) contains a pointer to the current Frame. Actions that want to alter the current Frame beyond just manipulating coordinates for subsequent Actions (e.g. changing the Frame size or adding velocity info etc) can do so via the **frm** variable (see e.g. Action\_Closest).

There are several variables needed for calculating the distance. First, we have two Vec3 classes (Vec3.h, which is already included from other headers) to store the XYZ coordinates of the points:

  Vec3 a1, a2;

If we are performing non-orthorhombic imaging we need to store the matrices which perform conversion from Cartesian to fractional coordinates and vice versa (using Matrix\_3x3 classes, Matrix\_3x3.h). Note that these are called ’ucell’ and ’recip’ respectively throughout CPPTRAJ, as these were the names used for the analogous structures in PTRAJ.

  Matrix\_3x3 ucell, recip;

Finally, we need a double to store the actual result of the distance calculation:

  double Dist;

In the first part of the actual calculation, we calculate the centers of the coordinates in Mask1\_ and Mask2\_, either mass-weighted or not depending on useMass\_, using the appropriate functions from the Frame class (Frame.h):

  if (useMass\_) {

    a1 = frm.Frm().VCenterOfMass( Mask1\_ );

    a2 = frm.Frm().VCenterOfMass( Mask2\_ );

  } else {

    a1 = frm.Frm().VGeometricCenter( Mask1\_ );

    a2 = frm.Frm().VGeometricCenter( Mask2\_ );

  }

Note that here we are using the Frm() function, which returns a constant (i.e. non-modifiable) reference to the current Frame; if we wanted to actually manipulate the coordinates we would have to call ModifyFrm().

Next, we get the distance between the coordinates stored in a1 and a2. For non-orthorhombic imaging we first need to convert the current box coordinates (stored in the Frame class in double precision as 3 lengths and 3 angles) into the coordinate conversion matrices using the ToRecip() function of the Box class (Box.h). Then, depending on the type of imaging that needs to be performed we call the appropriate distance calculation routine (DIST2\_XXX, found in DistRoutines.h):

  switch ( image\_.ImageType() ) {

    case NONORTHO:

      frm.Frm().BoxCrd().ToRecip(ucell, recip);

      Dist = DIST2\_ImageNonOrtho(a1, a2, ucell, recip);

      break;

    case ORTHO:

      Dist = DIST2\_ImageOrtho(a1, a2, frm.Frm().BoxCrd());

      break;

    case NOIMAGE:

      Dist = DIST2\_NoImage(a1, a2);       break;

  }

  Dist = sqrt(Dist);

Last, we add the result to the DataSet and return Action::OK. Since DataSet is just an interface we pass in the address of Dist (&Dist) to let the underlying DataSet framework take care of the fact that it is a double.

  dist\_-\>Add(frameNum, &Dist);

  return Action::OK;

}

### Print() - Any post-processing

The Print() function is called once all input frames have been read in, and is used if there is anything that should be printed outside the normal DataFile/DataSet framework (e.g. hydrogen bond averages in the hbond action) or if there are any additional calculations that need to be performed (e.g. finishing up matrix calculations in the matrix action). In this example we’re only calculating a simple distance; the output is handled by the DataFile/DataSet framework, so we implement a blank Print() function in the header:

void Print() {}

Add the Action to the Command class
-----------------------------------

Now that the class implementation is complete, we need to let cpptraj know how to call it. This is currently done using a “static” Class, Command (Command.cpp), which is initialized by the Cpptraj class via Command::Init() when the program starts. Command::Init() makes use of the Command::AddCmd() function to add the Command, set its destination, and any associated keywords. The Command::AddCmd() function looks like:

void Command::AddCmd(DispatchObject<span>\*</span> oIn, Cmd::DestType dIn, int nKeys, ...)

where **oIn** is a pointer to the DispatchObject (Exec-, Action-, Analysis-, or Deprecated-derived class), **dIn** determines how the Command will be processed, **nKeys** is the number of keywords associated with the command, and the remaining arguments are the command keywords. For example:

Command::AddCmd( new Action\_Rmsd(), Cmd::ACT, 2, rms, rmsd );

Adds a new instance of the Action-derived class Action\_Rmsd, sets its destination as Action (Cmd::ACT), and sets 2 associated command keys, “rms” and “rmsd”.

To make navigation of Commands.cpp easier, you can search for ACTION (or ANALYSIS if adding an Analysis) to go where things need to be added. First add the class to Commands.cpp with the appropriate ’\#include’. Includes should be in alphabetical order within their given section.

\#include Action\_Dihedral.h

**\#include Action\_Distance.h**

\#include Action\_Hbond.h

Then add the command to Command::Init() using Command::AddCmd(), e.g.:

Command::AddCmd( new Action\_Dipole(), Cmd::ACT, 1, dipole );

**Command::AddCmd( new Action\_Distance(), Cmd::ACT, 1, distance );**

Command::AddCmd( new Action\_DistRmsd(), Cmd::ACT, 2, drms, drmsd ); 
