function ax = lax(varargin)

%%
subrows = varargin{1};
subcols = varargin{2};

%%
len = length(varargin)-3;
sub = 1;
for i = 1 : len
    item = varargin{i+2};
    var = whos('item');
    if strcmp(var.class, 'char')
        title(ax(sub-1), item);
    else
        ax(sub) = subplot(subrows, subcols, sub);
        imshow(item, []);
        sub = sub + 1;
    end
end

if varargin{end}
    linkaxes(ax);
end