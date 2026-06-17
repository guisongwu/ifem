function y = mp(x)
% wrapper for missing multiple precison tool

if( isnumeric(x) )
    y = x;
elseif( isa(x,'char') )
    y = str2num(x);
else
    error('Unkown mp variable!!!\n');
end