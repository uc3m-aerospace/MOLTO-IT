%--------------------------------------------------------------------------0
function  molto_it_json(jsonfile)
%--------------------------------------------------------------------------
%	MOLTO-IT Software Computation Core										
%																			
%	This program is developed at the Universidad Carlos III de Madrid,		
%   as part of a PhD program.										
%																			
%   The software and its components are developed by David Morante González															
%																		
%   The program is released under the MIT License
%
%   Copyright (c) 2019 David Morante González															
%																			
%--------------------------------------------------------------------------
%
%    Function that call molto-it given a problem defined in a json formated
%    file
%
%--------------------------------------------------------------------------

input  = jsondecode(jsonfile);
[~] = molto_it(input);