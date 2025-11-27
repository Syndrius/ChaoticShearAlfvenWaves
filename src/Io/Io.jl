"""

Module for the reading and writing files using JLD2. Includes
 - Reading and writing our input data structures.
 - Writing solutions to file.

"""
module Io


using JLD2
using Printf
using EllipsisNotation

#good
include("FromFile.jl")

export inputs_from_file
export evals_from_file
export efunc_from_file
export surfaces_from_file

#good
include("ToFile.jl")

export inputs_to_file
export evals_to_file
export efuncs_to_file
export surfaces_to_file

end
