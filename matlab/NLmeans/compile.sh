#!/usr/bin/env  bash 
## 
#   Compile mex files for NLmeans
#
# - Michael Eager (michael.eager@monash.edu.au)
# - Monash Biomedical Imaging 
#
#  "$Id$"
#  "$Date$"
#
# Copyright (C) 2014 Michael Eager  (michael.eager@monash.edu)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
###############################################



if [ -x mex ]; then
cd $(dirname $0)/MRIDenoisingModified/
mex -v -compatibleArrayDims COMPFLAGS='$COMPFLAGS -O3'  OPTIM='-O3 -DNDEBUG' CXXOPTIMFLAGS='-O3 -DNDEBUG' LDOPTIMFLAGS='-O3 -DNDEBUG'  myODCT3d.cpp
mex -v -compatibleArrayDims COMPFLAGS='$COMPFLAGS -O3' OPTIM='-O3 -DNDEBUG'   CXXOPTIMFLAGS='-O3 -DNDEBUG' LDOPTIMFLAGS='-O3 -DNDEBUG' myMBONLM3D.cpp
mex -v -compatibleArrayDims COMPFLAGS='$COMPFLAGS -O3' OPTIM='-O3 -DNDEBUG'  CXXOPTIMFLAGS='-O3 -DNDEBUG' LDOPTIMFLAGS='-O3 -DNDEBUG' myRINLM3d.cpp
#mex -v -compatibleArrayDims COMPFLAGS='$COMPFLAGS -Ofast'  OPTIM='-Ofast -DNDEBUG' CXXOPTIMFLAGS='-Ofast -DNDEBUG' LDOPTIMFLAGS='-Ofast -DNDEBUG'  myODCT3d.cpp
#mex -v -compatibleArrayDims COMPFLAGS='$COMPFLAGS -Ofast' OPTIM='-Ofast -DNDEBUG'   CXXOPTIMFLAGS='-Ofast -DNDEBUG' LDOPTIMFLAGS='-Ofast -DNDEBUG' myMBONLM3D.cpp
#mex -v -compatibleArrayDims COMPFLAGS='$COMPFLAGS -Ofast' OPTIM='-Ofast -DNDEBUG'  CXXOPTIMFLAGS='-Ofast -DNDEBUG' LDOPTIMFLAGS='-Ofast -DNDEBUG' myRINLM3d.cpp
else

cd $(dirname $0)
matlab -nosplash -nodesktop -r "cd `pwd`/MRIDenoisingModified; \
mex -v -compatibleArrayDims COMPFLAGS='$COMPFLAGS -O3'  OPTIM='-O3 -DNDEBUG' CXXOPTIMFLAGS='-O3 -DNDEBUG' LDOPTIMFLAGS='-O3 -DNDEBUG'  myODCT3d.cpp;  \
mex -v -compatibleArrayDims COMPFLAGS='$COMPFLAGS -O3' OPTIM='-O3 -DNDEBUG'   CXXOPTIMFLAGS='-O3 -DNDEBUG' LDOPTIMFLAGS='-O3 -DNDEBUG' myMBONLM3D.cpp; \
mex -v -compatibleArrayDims COMPFLAGS='$COMPFLAGS -O3' OPTIM='-O3 -DNDEBUG'  CXXOPTIMFLAGS='-O3 -DNDEBUG' LDOPTIMFLAGS='-O3 -DNDEBUG' myRINLM3d.cpp; quit"

fi


