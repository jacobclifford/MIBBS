et=load('energyth.txt');
% the sampling distribution of the logodds is a normal distribution.  The expect mean is zero since the odds 
% should be one, hence logodds is zero, and std of the  sampling distribution is the standard error,
% which can be approximately calculated as sqrt of the sum
% of inverses of table counts.  
% the first column of 'a' contains the zscore for the spacer cutoff defined
% by the corresponding row in the third column of 'a'.  (the third column
% of 'a' is the spacer, second column is se, which isn't needed, since
% zscore is = logodds/se...  hence to calculate the logodds =
% a(:,1)*a(:,2)
% we want a one-sided pvalue: the probability of a table with dc having
% larger odds than the one attained by the simulation.
% since normcdf(x) returns the fraction of tables with logodds ratio
% smaller than the one attained, we obtain the desired value by subtracting
% that from one.
% hence the variable 'aa' below is the fraction of logodds ratios that are
% expected with a conditional probability of adjacent to Twist that is
% greater than that attained.
% since this may be a small number we plot the -log10 of this pvalue.
%aa = 1-normcdf(a(:,1));
%ee = [aa, a(:,3) ]
%ee = -log10( aa );
% len = max(a(:,2));
% plot(ee)
% xlabel('spacer cutoff in bp','FontSize',12,'FontName','Times')
% ylabel('-log10(pvalue of logodds ratio)','FontSize',12,'FontName','Times')
% legend('pwm threshold of pvalue = .0001','FontSize',12,'FontName','Times')
% print(gcf,'-depsc','logoddpred.eps')
 hold on
% i=1;
% j==0 is 10^-6
%j==1 is 10^-5
%j==2 is 10^-4
%j==3 is 10^-3
 for j = 0:4   %1:length(a(:,1))/len
     int2str(j)
     a =load(['da',int2str(j),'.txt']);
     aa = 1-normcdf(a(:,1));
%ee = [aa, a(:,3) ]
ee = -log10( aa );
%hold on
    %plot(ee([i:j*len]))
      if(j==0) plot(ee,'r') 
      end
      if(j==1) plot(ee,'y') 
      end
      if(j==2) plot(ee,'g')
          ee
      end
      if(j==3) plot(ee,'m') 
      end
      if(j==4) plot(ee,'k') 
      end
%     scatter(a([i:j*len],3),a([i:j*len],1))
%    % scatter(a([i:j*len],3),a([i:j*len],2),'r')
%     i =j*len;
%    % legend('CACATGT','SE','ACACAAA','SE')
 end
  %legend(['PWM pvalue=10^{-6}, fitnes thres=',int2str(et(1,:))], ['PWM pvalue=10^-5, fitnes thres=',int2str(et(2,:))],['PWM pvalue=10^-4, fitnes thres=',int2str(et(3,:))],['PWM pval=10^-3, fitnes thres=',int2str(et(1,:))])
   h=gca;
 set(h, 'FontSize', 18);  
  legend(['PWM pvalue=10^{-6}, fitness thres D_c=',num2str(et(1,2),2),' D_u='  num2str(et(1,3),2)], ['PWM pvalue=10^{-5}, fitness thres D_c=',num2str(et(2,2),2),' D_u=' num2str(et(2,3),2)],['PWM pvalue=10^{-4}, fitness thres D_c=',num2str(et(3,2),2),' D_u=' num2str(et(3,3),2)],...
        ['PWM pval=10^{-3}, fitness thres D_c=',num2str(et(4,2),2),' D_u=' num2str(et(4,3),2)],['PWM pval=10^{-0}, fitness thres D_c=',num2str(et(5,2),2),' D_u=' num2str(et(5,3),2)],'Location','NorthOutside')
  %  bb = xlabel({'fitness =$$ \epsilon$$ '});
   XTick = [1 2 3] ;
  XTickLabel ={'[0,30]bp','[31,60]bp','[61,90]bp'};
  set(h,'XTick',XTick,'XTickLabel',XTickLabel)
%set(bb,'interpreter','latex', 'FontSize', 12);
 legend('boxoff')
  xlabel('spacer windows','FontSize',18,'FontName','Times')
 ylabel('-log10(pvalue of logodds ratio)','FontSize',18,'FontName','Times')
 print(gcf,'-depsc','logoddsPredpvaAlltw.eps')
