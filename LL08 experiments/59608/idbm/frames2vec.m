function [ vec ] = frames2vec( frames, Ns, direction, window, synthesis )
% FRAMES2VEC Converts frames to signal using weighted overlap-and-add synthesis. 
%
%   A=FRAMES2VEC(B,N,D,W,S) converts frames in matrix B (stored as rows 
%   or columns as specified by D) into signal using weighted overlap-and-add
%   procedure specified by synthesis window W and synthesis type S.
%
%   Summary
%
%           B is a matrix of frames
%
%           N is a frame shift (in samples) originally used to split 
%             the signal into frames; alternatively, N can specify the matrix 
%             of indexes originally used to split the signal into frames
%             (i.e., the second output argument of the VEC2FRAMES function)
%
%           D specifies if the frames in B are rows or columns,
%             i.e., D = 'rows' or 'cols', respectively
%
%           W is the synthesis window function to be applied to each frame,
%             given as either a function handle, e.g., W = @hanning
%             or as a vector of window samples, e.g., W = hanning( M )
%
%           S specifies overlap-and-add synthesis of one of three types:  
%                   S = 'A&R';      % Allen & Rabiner's method
%                   S = 'G&L';      % Griffin & Lim's method
%                   S = 'Vanilla';  % Vanilla approach (no synthesis window)
%
%           A is the synthesized output signal 
%
%   Examples
%
%           % generate input signal samples
%           signal = [ 1:20 ]
%
%           % divide the input signal into seven-sample-long frames with a shift
%           % of three samples, pad the last frame with zeros so that no samples 
%           % are discarded, apply the Hanning analysis window to each frame and
%           % return frames as columns vectors
%           frames = vec2frames( signal, 7, 3, 'cols', @hanning, 0 )
%
%           % synthesize the frames back to signal using weighted overlap-and-add
%           % procedure; note the extra padding on the reconstructed signal
%           reconstructed = frames2vec( frames, 3, 'cols', @hanning, 'G&L' )
% 
%   See also VEC2FRAMES.

%   Author: Kamil Wojcicki, UTD, July 2011


    % usage information
    usage = 'usage: [ vec ] = frames2vec( frames, indexes, direction, window, synthesis );';

    % default settings 
    switch( nargin )
    case { 0, 1, 2 }, error( usage );
    case 3, synthesis='G&L'; window=@hanning; 
    case 4, synthesis='G&L'; 
    end

    % input validation
    if( isempty(frames) || isempty(Ns) || isempty(direction) ), error( usage ); end;
    if( Ns<1 ), error( usage ); end;


    % process based on frame direction
    switch( direction )

    % rows as frames
    case 'rows'     

        % get number of frames and frame length
        [ M, Nw ] = size( frames );

        if( length(Ns)==1 )
            % if frame duration provided, generate framing indexes
            indf = Ns*[ 0:(M-1) ].';                          % indexes for frames      
            inds = [ 1:Nw ];                                  % indexes for samples
            indexes = indf(:,ones(1,Nw)) + inds(ones(M,1),:); % framing indexes
        else 
            % otherwise, framing indexes were provided
            indexes = Ns;
            Ns = indexes(2,1)-indexes(1,1);
        end

        % determine signal duration
        L = max(indexes(:));

        % if synthesis window provided as a function handle
        % then generate synthesis window samples
        if( isa(window,'function_handle') )
            window = window( Nw );
        end
        window = window(:).';

        % allocate storage
        vec = zeros(1, L); 
        wsum = zeros(1, L); 

        % overlap-and-add syntheses 
        switch(upper(synthesis))
    
        % Allen & Rabiner's method
        case {'ALLEN & RABINER','A&R'}                   
    
            % overlap-and-add frames
            for m = 1:M, vec(indexes(m,:)) = vec(indexes(m,:)) + frames(m,:); end;

            % overlap-and-add window samples
            for m = 1:M, wsum(indexes(m,:)) = wsum(indexes(m,:)) + window; end;

            % for some tapered analysis windows, endpoint samples are very close 
            % to zero; as a consequence of this, wsum can be very close to
            % zero at its endpoints; division of vec by wsum at points where 
            % wsum is close to zero can produce large impulses; if you are 
            % experiencing this issue, one approach to address this is to limit
            % the lower bound of wsum as follows:
            %
            % wsum( wsum<1E-2 ) = 1E-2;  

            % divide out summed-up analysis windows
            vec = vec./wsum;
    
        % Griffin & Lim's method
        case {'GRIFFIN & LIM','G&L'}
    
            % apply synthesis window
            frames = frames * diag( window );

            % overlap-and-add frames
            for m = 1:M, vec(indexes(m,:)) = vec(indexes(m,:)) + frames(m,:); end;

            % overlap-and-add squared window samples
            for m = 1:M, wsum(indexes(m,:)) = wsum(indexes(m,:)) + window.^2; end;

            % for some tapered analysis windows, endpoint samples are very close 
            % to zero; as a consequence of this, wsum can be very close to
            % zero at its endpoints; division of vec by wsum at points where 
            % wsum is close to zero can produce large impulses; if you are 
            % experiencing this issue, one approach to address this is to limit
            % the lower bound of wsum as follows:
            %
            % wsum( wsum<1E-2 ) = 1E-2;  

            % divide out squared and summed-up analysis/synthesis windows
            vec = vec./wsum;

        % vanilla approach
        case {'VANILLA'}

            % overlap-and-add frames
            for m = 1:M, vec(indexes(m,:)) = vec(indexes(m,:)) + frames(m,:); end;
    
        % unsupported approach
        otherwise

            error(sprintf('%s: synthesis type not supported.', synthesis));
    
        end 


    % columns as frames
    case 'cols'

        % get frame length and the number of frames 
        [ Nw, M ] = size( frames );

        if( length(Ns)==1 )
            % if frame duration provided, generate framing indexes
            indf = Ns*[ 0:(M-1) ];                            % indexes for frames      
            inds = [ 1:Nw ].';                                % indexes for samples
            indexes = indf(ones(Nw,1),:) + inds(:,ones(1,M)); % framing indexes
        else
            % otherwise, framing indexes were provided
            indexes = Ns;
            Ns = indexes(1,2)-indexes(1,1);
        end

        % determine signal duration
        L = max(indexes(:));

        % if synthesis window provided as a function handle
        % then generate synthesis window samples
        if( isa(window,'function_handle') )
            window = window( Nw );
        end
        window = window(:);

        % allocate storage
        vec = zeros(L, 1); 
        wsum = zeros(L, 1); 
       
        % overlap-and-add syntheses 
        switch(upper(synthesis))
    
        % Allen & Rabiner's method
        case {'ALLEN & RABINER','A&R'}                   
    
            % overlap-and-add frames
            for m = 1:M, vec(indexes(:,m)) = vec(indexes(:,m)) + frames(:,m); end;

            % overlap-and-add window samples
            for m = 1:M, wsum(indexes(:,m)) = wsum(indexes(:,m)) + window; end;

            % for some tapered analysis windows, endpoint samples are very close 
            % to zero; as a consequence of this, wsum can be very close to
            % zero at its endpoints; division of vec by wsum at points where 
            % wsum is close to zero can produce large impulses; if you are 
            % experiencing this issue, one approach to address this is to limit
            % the lower bound of wsum as follows:
            %
            % wsum( wsum<1E-2 ) = 1E-2;  

            % divide out summed-up analysis windows
            vec = vec./wsum;
    
        % Griffin & Lim's method
        case {'GRIFFIN & LIM','G&L'}
    
            % apply synthesis window
            frames = diag( window ) * frames;

            % overlap-and-add frames
            for m = 1:M, vec(indexes(:,m)) = vec(indexes(:,m)) + frames(:,m); end;

            % overlap-and-add squared window samples
            for m = 1:M, wsum(indexes(:,m)) = wsum(indexes(:,m)) + window.^2; end;

            % for some tapered analysis windows, endpoint samples are very close 
            % to zero; as a consequence of this, wsum can be very close to
            % zero at its endpoints; division of vec by wsum at points where 
            % wsum is close to zero can produce large impulses; if you are 
            % experiencing this issue, one approach to address this is to limit
            % the lower bound of wsum as follows:
            %
            % wsum( wsum<1E-2 ) = 1E-2;  

            % divide out squared and summed-up analysis/synthesis windows
            vec = vec./wsum;
    
        % vanilla approach
        case {'VANILLA'}

            % overlap-and-add frames
            for m = 1:M, vec(indexes(:,m)) = vec(indexes(:,m)) + frames(:,m); end;
    
        % unsupported approach
        otherwise

            error(sprintf('%s: synthesis type not supported.', synthesis));

        end

    end


% EOF 
