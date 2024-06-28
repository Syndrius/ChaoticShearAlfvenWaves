"""

Module for the reading and writing files. Includes
 - Reading and writing our input data structures.
 - Writing solutions to file.

"""
module Io

using MID.Inputs
using MID.Geometry

include("FromFile.jl")

export inputs_from_file


include("ToFile.jl")

export inputs_to_file
export eigfuncs_to_file
export eigvals_to_file


end