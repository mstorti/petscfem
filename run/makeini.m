flag=xnod(:,1)>0;

ulh= 1;
ulu=.4;

urh = .8;
uru = .4;

u = (1-flag) * [ulh*ulu 0 ulh] + flag * [urh*uru 0 urh];
asave("channel_step.ini",u);
