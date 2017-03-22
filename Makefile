#
# Metagenomic Assignment
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
# written by Conrad Shyu (conradshyu@hotmail.com)
#
# Author's comments:
# ------------------
# Center of the Study of Biological Complexity (CSBC)
# Department of Microbiology and Immunology
# Virginia Commonwealth University
# Richmond, VA 23298
#
# revised on March 18, 2014
#
all: samfile assign

samfile:
	g++ -I. -O3 samfile.cpp -o samfile -fopenmp

assign:
	g++ -I. -O3 assign.cpp strain.cpp species.cpp -o assign -fopenmp

clean:
	rm -f samfile assign
