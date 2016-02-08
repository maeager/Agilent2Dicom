#!/usr/bin/env  bash 
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


source ~/.bashrc; 
if test ${MASSIVE_USERNAME+defined}; 
then 


module purge
module load build gmp mpc mpfr
module virtualgl matlab git
fi

cd $(dirname $0)
git clone https://github.com/vlfeat/vlfeat.git
cd vlfeat
make
matlab -nodesktop -nosplash -r 'cd ${PWD}/toolbox;run  vl_setup.m; quit'
