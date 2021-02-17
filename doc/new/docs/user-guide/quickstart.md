# Quick Start

## Workflow

### Starting from a TGMM config file

The commands listed below refer to executables that are run from a terminal
shell (`cmd` or `powershell` on Windows, `sh` on Linux).

#### 1. Compute a hierarchical segmentation:

        ProcessStack <configFile> <frame>

Computes the hierarchical segmentation, truncating the tree at ``minTau``.
Produces two files

Relevant parameters:


Parameter          | Description
------------------ | -------------------------------------------------------
imgBasename        | Path to a tiff or klb with the stack for a single timepoint.
radiusMedianFilter | The full width of the filter will be ``2*radiusMedianFilter+1``
minTau             | The agglomeration threshold at which to cut the tree. See TGMM_UserGuide.pdf for details.
backgroundThr      | Supervoxels are only considered for foreground regions.
conn3D             | Either 6,8,74.  Specifies how many neighboring voxels are considered connected.  74 is a 5x5x3 region of neighbors.
   


#### 2. (optional) Output a segmentation:

        ProcessStack <binFile> <tau> <minSuperVoxelSzPx>

Will generate a segmentation and output the result to a klb.

#### 3. Run Tracking:

        TGMM <configFile> <startTime <endTime>

### The config file

[**Example**](example-config.md)

[**Parameters Guide**](parameters.md)

All parameters, as well as the location of the input and output data is stored
in a single configuration file. 

There are three sections:

* Path parameters
* Main parameters
* Advanced parameters

Default values for any missing parameters.

An "experiment log" text file is saved to the output folder.  It records all
the settings used including the path settings. 

## Tips

### ProcessStack

* See the tgmm paper and the TGMM_UserGuide.pdf for detailed descriptions of 
  everything.
* Avoid long path names (full path must be less than 260 characters).
* ProcessStack outputs data in the same folder as the input data, but using a
  slightly modified name.
* ProcessStack often wants input filenames with the file extension removed.
  So instead of ``path/to/my_timepoint.tif`, use `path/to/my_timepoint`.
* When specifying paths on windows, it's fine to use (unix-style) 
  backslashes: `C:/path/to/my/files`. Windows normally uses forward slashes.
  On modern versions of windows, backslashes in paths get translated
  appropriately.  If you really want to use forward slashes in the config
  file, they need to be escaped with another slash: `C:\\path\\to\\my\\files`.
  Either approach should work. 

    