function mat2scr(x, fmt, label)

% Function mat2scr(x, fmt, label)
% Function to display a matrix as rows of ascii text.
% Similar to mat2str but sends result directly to the screen.
% fmt controls sprintf - ' %.0f' gives integers with single spaces.
% If fmt is not given, '%6.2f' is used by default.
% If x is complex, the imag parts are printed on rows below each row
% of real parts, and annotated '+j'.
% If label is given, this text is used as a label.

if isempty(x), return, end

if nargin < 2,
  fmt = '%6.2f';
end

if nargin > 2,
  disp(label);
end

[m,n] = size(x);

for j = 1:m,
  disp(sprintf(fmt, real(x(j,:))));
  if ~isreal(x),
    disp(['+j' sprintf(fmt, imag(x(j,:)))]);
  end
end

return
