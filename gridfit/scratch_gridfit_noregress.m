% scratch_gridfit_noregress.m

all_gr = linspace(0,2,25);
all_gc = linspace(0.1,100,25);
all_b  = linspace(-0.5,0.5,25);
all_n  = linspace(0.5,10,25);

[gr1, gr2 gc1 gc2 b1 b2 n] = ndgrid(all_gr,all_gr,all_gc,all_gc,all_b,all_b,all_n);

all_params = [reshape(gr1,numel(gr1),1) reshape(gr2,numel(gr1),1) reshape(gc1,numel(gr1),1) reshape(gc2,numel(gr1),1) reshape(b1,numel(gr1),1) reshape(b2,numel(gr1),1) reshape(n,numel(gr1),1)];
