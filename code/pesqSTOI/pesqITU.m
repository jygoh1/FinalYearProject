function pmos = pesqITU(fs,ref_file,deg_file,swap)
% PESQ - Perceptual evaluation of speech quality. 
%   
%   Useage:     
%       pmos = pesq(fs,ref_file,deg_file,swap);
%
%   Inputs:     
%       fs      sampling frequency (8 or 16 kHz)
%       ref     reference speech file   (wav, 16-bits)
%       deg     degraded speech file    (wav, 16-bits)
%       swap    swap byte order (TRUE 1 or FALSE 0, default = 0)
%
%   Outputs:
%       pmos    the predicted opinion score from pesq
%

%
%   Notes: 
%       1. PESQ is the ITU-T standard recommendation P.862 for objective perceptual 
%       evaluation of narrow-band tlephone speech quality. It is mostly applied
%       for assessment of narrow-band telephone networks and speech codecs.
%   
%       2. MATLAB MEX-version using the ITU-T C source
%           pesq_itu.mexmac     compiled and tested on OS X 10.4
%           pesq_itu.dll
%
%   References:
%   
%   
%   TODO:
%
%**************************************************************************
%	Copyright (C) CLEAR 2008
%   Version: 
%   Author:     Nick Gaubitch (MEX-version using ITU-T C source)
%   Date:       Mar 2008
%**************************************************************************
if (nargin < 4)  swap = 0; end
if (fs ~= 8000 & fs ~= 16000) error('Invalid sampling frequncy. The sampling frequency must be 8 or 16kHz'); end
pmos = pesq_itu(fs, ref_file, deg_file, swap);