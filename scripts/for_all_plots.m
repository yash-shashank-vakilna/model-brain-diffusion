names=dir('*.csv');
for i=1:length(names)
    format_hs_graph(names(i).name);
end
