 a=load('energy.txt');
 a2=load('energy2.txt');
 dc = load('dcondit.txt');
 du = load('ducondit.txt');
%aaa=load('t.txt');
 c=ranksum(dc,du);
%  [g,c]=ttest2(dc,du);
bin = zeros(1,300);
bin2 = zeros(1,300);
[aa,bb ] =size(a(:,1));

bin3 = zeros(1,aa);

binmax = 20;
%for i = 1: size(a(:,1))
 for i = 1:size(a(:,1))
  
x = a2( i  , : ) ;
y =  a( i , : ) ;
 b =  ranksum( y, x) ;
% [gg, b] =  ttest2( y, x) ;
 bin3(i) = b;
  binum = -10*log10(b) ;
%binum = b*10000 ;
 binum = ceil( binum) +1 ;
 bin(binum) = bin(binum) + 1;
  if(binum > 30 ) 
      bin2(binum) = bin2(binum) + 1;
 end
 end

bar(bin,'b')
 binum = -log10(c)*10 ;
  binum = ceil( binum) +1 ;
  sum(bin)
  hold on

  binww = zeros(1,300);
  
 binww(binum) = binww(binum) +50;
bar(binww,'r')  % this was 'r' but modified since color red wasn't displaying for small figures

%legend('permutations','DC vs DU');
%hold off
%bar(bin2)
%h = findobj(gca,'Type','patch');
%set(h,'FaceColor','r','EdgeColor','w')
  
%set(get(gca,'Y'),'interpreter','latex', 'FontSize', 25);
 %   set( h ,'interpreter','latex', 'FontSize', 8);
 %f=plot(A);
 % str2=['saveas(h, ''du' int2str(aaa) '.eps'',''eps'')'];
 % eval(str2);
%hold off
set(gca,'FontSize',18)
%set(gca,'xscale','log');
%title('Permutation test for checking significance between subpopulations of Dorsal sites','FontSize',12)
xlabel('-10 * log_{10}(pvalue of ranksum test) ')
ylabel('frequency of pvalue ','FontSize',18)
 print(gcf,'-depsc','rs.eps')