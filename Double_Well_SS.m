%*********************************************SS Double Well VNA ER Updated***************************************************
%test_flag=1, do ER only; flag=2, do VNA only; flag=3, do VNA then ER; flag=4, do ER then VNA

function [Data1, AvgPwr1, e1, e1_pat, Data2, AvgPwr2, e2, e2_pat]= Double_Well_SS(FFL_or_Phip, phin_path1, phin_sweep1, phip_levels, power, avg, bw, fstart, fstop, npoints, sample_range, pattest, ntest, edges, bitwise, test_flag, tag, plot_flag)
    global imb iob icb itb idb ibs ffl
    global PATTERN pat1 pat1b pat2b pat3b pat4b pat100b pat400b pat1000b pat4000b pat10000b pat40000b pcl
    global imb_chan iob_chan itb_chan idb_chan icb_chan ibs_chan ffl_chan
    global v gnd fridge chassis
    %Sweep parameers
    if nargin<1 || isempty(FFL_or_Phip),  FFL_or_Phip=  "Phip"; end
    if nargin<2 || isempty(phin_path1),  phin_path1=  0:.01:0.4; end
    if nargin<3 || isempty(phin_sweep1), phin_sweep1= 0.4:-.005:-0.4; end
    if nargin<4 || isempty(phip_levels), phip_levels=  [0,-0.4]; end
    %VNA parameters
    if nargin<5 || isempty(power), power= -45; end
    if nargin<6 || isempty(avg), avg= 2; end
    if nargin<7 || isempty(bw), bw=1e3; end
    if nargin<8 || isempty(fstart), fstart= 6.7e9; end
    if nargin<9 || isempty(fstop), fstop=6.7000001e9; end
    if nargin<10 || isempty(npoints), npoints= 31; end
    if nargin<11 || isempty(sample_range), sample_range=[fstart, fstop]; end
    %Comp parameters
    if nargin<12 || isempty(pattest), pattest=1; end
    if nargin<13 || isempty(ntest), ntest=1; end    %ntest is number of consecutive single shot tests to perform in each state
    if nargin<14 || isempty(edges), edges=1; end
    if nargin<15 || isempty(bitwise), bitwise=true; end
    %Label/Flag
    if nargin<16 || isempty(test_flag), test_flag=2; end
    if nargin<17 || isempty(tag), tag=""; end
    if nargin<18 || isempty(plot_flag), plot_flag=1; end
    %Global equipment definitions/ error catches
    if test_flag~=1 && (~exist('chassis', 'var') || isempty(chassis))
        chassis=input("Enter chassis number of Keysight VNA (usually 3 for ADRs): ");
    end
    if ~exist('fridge', 'var') || isempty(fridge)
        fridge=input("Enter fridge in use, for temp grabs ('adr1', 'adr2', 'He', 'dr', etc): ");
    end
    if any(sample_range < 1e6, 'all') || any([fstart, fstop] < 1e6, 'all')
        disp("Frequencies may not be defined correctly, less than 1 MHz");
    end

    %Comp initializing
    e1=[];
    e2=[];
    e1_pat=[];
    e2_pat=[];
    e1_ns=[]; 
    e1_ds=[];
    e1_os=[];
    e2_ns=[];
    e2_ds=[];
    e2_os=[];
    %splitcurrent = 0; %May be needed in future edits

    if test_flag~=2
        bits=sum(PATTERN(pattest).vdat,"all");
        if isnumeric(bitwise) && bitwise>0
            bits=round(bitwise);
        end
    end

    %VNA initializing
    f1={};
    Data1={};
    AvgPwr1={};
    f2={};
    Data2={};
    AvgPwr2={};   
    
    if test_flag ~= 1
        clear v
        v = PXI_VNA();
        v.setOutput(1);
        v.setSparam("S21");
        v.setPower(power)
        v.setAverage(avg)
        v.setBandwidth(bw)
        v.setStopFrequency(40e9)
        v.setStartFrequency(fstart)
        v.setStopFrequency(fstop)
        v.setNpts(npoints)
        v.getOutput
        v.getSparam
        v.getPower
        v.getAverage
        v.getBandwidth
        v.getStartFrequency
        v.getStopFrequency
        v.getNpts
        
        rang=linspace(fstart,fstop,npoints);  %Live results for first sample range only (if there are more than one)
        mini=sample_range(1,1);
        maxi=sample_range(1,2);
        idx1=find(rang>=mini,1,"first");
        idx2=find(rang<=maxi,1,"last");
        if isempty(idx1), idx1=1; end
        if isempty(idx2), idx2=length(rang); end
    end

    % icb is array of [DCFlux1, DCFlux2], also called icb(1) and icb(2), they are currents applied to induce flux in JJ1 and JJ2. 
    % Phi+= DCFlux1 - DCFlux2, Phi-= DCFlux1 + DCFlux2
    % DCFlux1 = (Phi- + Phi+)/2, DCFlux2 = (Phi- - Phi+)/2

    phin0=phin_path1(1);                    %Starting Phi minus for both sweeps
    phin_path2=2*phin_path1(1)-phin_path1;  %Define return initialization sweep based on first initialization sweep
    phin_sweep2=2*phin_sweep1(end)-phin_sweep1;%Define return test sweep based on first test sweep

    if length(icb_chan)==1                  %With 1 icb_chan bias, assume 1-line Phi-. Also use FFl 
        Phin=1;
        FFL_or_Phip="FFL";
        Phin_string=string("1-line Phi-");
    elseif length(icb_chan)==2              %With 2 icb_chan bias, assume 2-line Phi- = DCFLux1 + DCFlux2. Can be used with Phi+= DCFLux1 - DCFlux2 or with ffl_chan
        Phin=0;
        Phin_string=string("2-line Phi-");
    end

    %Set test to use Phip or FFL; if using fast flux, redefine terms to use ffl0 and ffl_PI and set all phip variables to 0
    if istrue(FFL_or_Phip) || contains(FFL_or_Phip, "f", "IgnoreCase", true)
        FF=1;                            
        ffl0=phip_levels(1); 
        ffl_PI=phip_levels(2);
        phip0=0;
        phip_PI=0;
        PI_string=string(ffl_PI);
    elseif isfalse(FFL_or_Phip) || contains(FFL_or_Phip, "p", "IgnoreCase", true)
        FF=0;
        phip0=phip_levels(1);   %Phi plus corresponding to 0 state
        phip_PI=phip_levels(2); %Phi plus corresponding to PI state
        PI_string=string(phip_PI);
    else
        disp("Specify 'Phip' or 'FFL' as first variable")
        return
    end

    %Each single shot test in ntest runs through full paths, similar to loop_times
    disp("************************************************************Sweep1")
    [e1, e1_pat, e1_ns, e1_ds, e1_os, Data1,  f1]=DW_Test(phin_path1,phin_sweep1);    
    disp("************************************Sweep2")
    [e2, e2_pat, e2_ns, e2_ds, e2_os, Data2,  f2]=DW_Test(phin_path2,phin_sweep2); 
    
    %turn off Phi+, Phi-, and FFL
    if Phin, c_ma(icb_chan,0); else, c_ma(icb_chan,[0,0]); end
    if FF, c_ma(ffl_chan,0); end   
    
    disp("***********************************************************************")
    %temp grab
    if contains(fridge, "he", "IgnoreCase", true)
        temp_string="4.2K He ";
    elseif strcmpi(fridge,"adr1") || strcmpi(fridge,"adr")
        temp_string=string(round(adr_temp*1000))+"mK ";
    elseif strcmpi(fridge,"adr2")
        temp_string=string(round(adr2_temp*1000))+"mK ";
    else
        temp_string=" ";
    end

    %Add tag for DUT based in icb_chan name
    icb_char=char(string(icb_chan(1)));
    und_idx=find(icb_char == '_', 1);
    if ~isempty(und_idx), dev_string=string(icb_char(1:und_idx-1)); else, dev_string=""; end

    %Add tag to device string labeling test order
    switch test_flag
        case 1
            dev_string=string(dev_string)+ " ER ";
        case 2
            dev_string=string(dev_string)+ " VNA ";
        case 3
            dev_string=string(dev_string)+ " VNA,ER ";
        case 4
            dev_string=string(dev_string)+ " ER,VNA ";
    end

    if plot_flag
        %Plot DW of comp BER
        if test_flag~=2
            cutoff=.5;
            e1_cutoff=e1>cutoff;
            e2_cutoff=e2>cutoff;  
            if (ischar(bitwise) || isstring(bitwise)) && contains(lower(string(bitwise)),"all")
                bit_visibility="on";
                pattern_visibility="on";
                cutoff_visibility="on";
                ds_os_visibility="on";
                plot_type="All ERs";
            elseif (ischar(bitwise) || isstring(bitwise)) && contains(lower(string(bitwise)),"cutoff")
                bit_visibility="off";
                pattern_visibility="off";
                cutoff_visibility="on";
                ds_os_visibility="off";
                plot_type="Pat Cutoff ER";
            elseif ((isnumeric(bitwise) && bitwise<=0) || ...
                    (islogical(bitwise) && ~bitwise)) || ...
                    ((ischar(bitwise) || isstring(bitwise)) && contains(lower(string(bitwise)),"pat"))
                bit_visibility="off";
                pattern_visibility="on";
                cutoff_visibility="off";
                ds_os_visibility="off";
                plot_type="Pat ER";
            else
                bit_visibility="on";
                pattern_visibility="off";
                cutoff_visibility="off";
                ds_os_visibility="off";
                plot_type="Bit ER";
            end
            % Code for plot of DW comp Total bit ER and also ER due to double switch, offcycle switch               
            
            leg1="Right Well "+string(mean(e1)*100)+"% Bit ER";
            leg2="Left Well " +string(mean(e2)*100)+"% Bit ER";
            leg3="Right Well "+string(mean(e1_ds+e1_os)*100)+"% DS+OS BER";
            leg4="Left Well " +string(mean(e2_ds+e2_os)*100)+"% DS+OS BER";
            leg5="Right Well "+string(mean(e1_pat)*100)+"% Pat ER";
            leg6="Left Well " +string(mean(e2_pat)*100)+"% Pat ER";
            leg7="Right Well "+string(mean(e1_cutoff)*100)+"% Pat ER >"+string(100*cutoff)+"% Cutoff";
            leg8="Left Well " +string(mean(e2_cutoff)*100)+"% Pat ER >"+string(100*cutoff)+"% Cutoff";
            figure
            hold on
            plot([1:ntest],e1,"o","DisplayName", string(leg1),"Visible",bit_visibility,"HandleVisibility",bit_visibility);
            plot([1:ntest],e2,'*',"DisplayName", string(leg2),"Visible",bit_visibility,"HandleVisibility",bit_visibility);
            plot([1:ntest],e1_ds+e1_os,"o","DisplayName", string(leg3),"Visible",ds_os_visibility,"YLimInclude","off","HandleVisibility",ds_os_visibility) 
            plot([1:ntest],e2_ds+e2_os,"*","DisplayName", string(leg4),"Visible",ds_os_visibility,"YLimInclude","off","HandleVisibility",ds_os_visibility)
            plot([1:ntest],e1_pat,"o","DisplayName",string(leg5),"Visible",pattern_visibility,"YLimInclude","off","HandleVisibility",pattern_visibility)
            plot([1:ntest],e2_pat,'*',"DisplayName",string(leg6),"Visible",pattern_visibility,"YLimInclude","off","HandleVisibility",pattern_visibility)
            plot([1:ntest],e1_cutoff,"o","DisplayName",string(leg7),"Visible",cutoff_visibility,"YLimInclude","off","HandleVisibility",cutoff_visibility)
            plot([1:ntest],e2_cutoff,'*',"DisplayName",string(leg8),"Visible",cutoff_visibility,"YLimInclude","off","HandleVisibility",cutoff_visibility)

            ylim([-.1,1.1]);
            xlabel("ntest")
            ylabel("Error Rate")
            legend('Location','northeast','Orientation','vertical');
            fig_name="SS "+plot_type+" at "+string(FFL_or_Phip)+"="+PI_string+" Phi-=" + string(phin_sweep1(end))+ " mA, "+ temp_string+string(dev_string)+string(ntest)+"x"+string(bits)+"b";    
            title({string(fig_name),string(tag)})
            saveas(gcf,string(fig_name) + " " +string(tag)+".fig");
            hold off
            close(gcf)
            
            bins = linspace(0,1,max(101,bits+1));
            figure 
            hold on
            histogram(e1,bins,'DisplayStyle','stairs','EdgeColor','b','LineWidth',3,"DisplayName",string(leg1),"Visible",bit_visibility,"HandleVisibility",bit_visibility);
            histogram(e2,bins,'DisplayStyle','stairs','EdgeColor','r','LineWidth',3,"DisplayName",string(leg2),"Visible",bit_visibility,"HandleVisibility",bit_visibility);
            histogram(e1_ds+e1_os,bins,"DisplayName",string(leg3),"Visible",ds_os_visibility,"HandleVisibility",ds_os_visibility)
            histogram(e2_ds+e2_os,bins,"DisplayName",string(leg4),"Visible",ds_os_visibility,"HandleVisibility",ds_os_visibility)
            histogram(e1_pat,bins,"DisplayName",string(leg5),"Visible",pattern_visibility,"HandleVisibility",pattern_visibility)
            histogram(e2_pat,bins,"DisplayName",string(leg6),"Visible",pattern_visibility,"HandleVisibility",pattern_visibility)
            histogram(e1_cutoff,bins,"DisplayName",string(leg7),"Visible",cutoff_visibility,"HandleVisibility",cutoff_visibility)
            histogram(e2_cutoff,bins,"DisplayName",string(leg8),"Visible",cutoff_visibility,"HandleVisibility",cutoff_visibility)

            set(gca,'YScale','log')
            xlabel('Error Rate by Test')
            ylabel('Count (Log Scale)')
            legend('Location','north','Orientation','vertical');
            fig_name="Hist SS "+plot_type+", "+string(FFL_or_Phip)+"="+PI_string+" Phi-=" + string(phin_sweep1(end))+ " mA, "+ temp_string+string(dev_string)+string(ntest)+"x"+string(bits)+"b";    
            title({string(fig_name),string(tag)})
            xlim([-.05,1.05])
            saveas(gcf,string(fig_name) + " " +string(tag)+".fig");
            close(gcf)
        end %test_flag
    end %plot_flag

    %Plot DW of VNA S21 at specified frequency(ies)
    if test_flag~=1
        rang=linspace(fstart,fstop,npoints);
        for k=1:size(sample_range,1)            
            mini=sample_range(k,1);
            maxi=sample_range(k,2);
            idx1=find(rang>=mini,1,"first");
            idx2=find(rang<=maxi,1,"last");
            if isempty(idx1) || isnan(idx1), idx1=1; end
            if isempty(idx2) || isnan(idx2), idx2=length(rang); end
            AvgPwr1={};
            AvgPwr2={};
            for n=1:ntest
                AvgPwr1{n,k}=[mean(Data1{n}(idx1:idx2,2)),   mean(atan(Data1{n}(idx1:idx2,4)./Data1{n}(idx1:idx2,3))),    mean(Data1{n}(idx1:idx2,3)),    mean(Data1{n}(idx1:idx2,4))];
                AvgPwr2{n,k}=[mean(Data2{n}(idx1:idx2,2)),   mean(atan(Data2{n}(idx1:idx2,4)./Data2{n}(idx1:idx2,3))),    mean(Data2{n}(idx1:idx2,3)),    mean(Data2{n}(idx1:idx2,4))];
            end
        end

        if plot_flag
            figure
            hold on
            for k=1:size(sample_range,1)     
                p1=cell2mat(AvgPwr1(:,k));
                p2=cell2mat(AvgPwr2(:,k));
    
                leg1="Right Well "+string( mean(p1(:,1:4:end)))+" Avg dB @";
                leg2="Left Well "+ string( mean(p2(:,1:4:end)))+" Avg dB @";
                leg3="Right Well "+string( mean(p1(:,2:4:end)))+" Avg Phase @";
                leg4="Left Well "+ string( mean(p2(:,2:4:end)))+" Avg Phase @";
                leg5="Right Well "+string( mean(p1(:,3:4:end)))+" Avg Re(S21) @";
                leg6="Left Well "+ string( mean(p2(:,3:4:end)))+" Avg Re(S21) @";
                leg7="Right Well "+string( mean(p1(:,4:4:end)))+" Avg Im(S21) @";
                leg8="Left Well "+ string( mean(p2(:,4:4:end)))+" Avg Im(S21) @";
                   
                plot([1:ntest], p1(:,1:4:end),"o","DisplayName", string(leg1)+string(round(mean(sample_range(k,:))./1e9,2)) + "GHz")
                plot([1:ntest], p2(:,1:4:end),'*',"DisplayName", string(leg2)+string(round(mean(sample_range(k,:))./1e9,2)) + "GHz")
                plot([1:ntest], p1(:,2:4:end),"o","DisplayName", string(leg3)+string(round(mean(sample_range(k,:))./1e9,2)) + "GHz", "Visible", "off", "YLimInclude","off","HandleVisibility","off")       
                plot([1:ntest], p2(:,2:4:end),'*',"DisplayName", string(leg4)+string(round(mean(sample_range(k,:))./1e9,2)) + "GHz", "Visible", "off", "YLimInclude","off","HandleVisibility","off") 
                plot([1:ntest], p1(:,3:4:end),"o","DisplayName", string(leg5)+string(round(mean(sample_range(k,:))./1e9,2)) + "GHz", "Visible", "off", "YLimInclude","off","HandleVisibility","off") 
                plot([1:ntest], p2(:,3:4:end),'*',"DisplayName", string(leg6)+string(round(mean(sample_range(k,:))./1e9,2)) + "GHz", "Visible", "off", "YLimInclude","off","HandleVisibility","off") 
                plot([1:ntest], p1(:,4:4:end),"o","DisplayName", string(leg7)+string(round(mean(sample_range(k,:))./1e9,2)) + "GHz", "Visible", "off", "YLimInclude","off","HandleVisibility","off") 
                plot([1:ntest], p2(:,4:4:end),'*',"DisplayName", string(leg8)+string(round(mean(sample_range(k,:))./1e9,2)) + "GHz", "Visible", "off", "YLimInclude","off","HandleVisibility","off")                 
            end
            axis tight
            hold off
            xlabel("ntest")
            ylabel("S21")
            legend
            fig_name="SS S21 dB at "+string(FFL_or_Phip)+"="+PI_string+" Phi-=" + string(phin_sweep1(end))+ " mA, "+temp_string+string(dev_string)+string(round(mean(sample_range(1,:))./1e9,2)) + "GHz";
            title({string(fig_name),string(tag)})
            saveas(gcf,string(fig_name) + " " +string(tag)+".fig")
            close(gcf)
    
            %If frequency span more than 1MHz, plot Heatmaps of DW S21 vs Phi- & Freq
            if abs(fstop-fstart) >.001e9      
                div=ceil((npoints-1)/100)*10+1;
                ylabels=string(round(f1{1}'./1e6));
                ylabels(  mod(1:length(f1{1}),div)~=1   ) ="";
                xlabels1=string([1:ntest]);
                xlabels1(  mod(1:ntest,10)~=0   ) ="";
                c1=cell2mat(Data1(:)');
                xlabels2=string([1:ntest]);
                xlabels2(  mod(1:ntest,10)~=0   ) ="";
                c2=cell2mat(Data2(:)');
            
                figure
                h=heatmap([1:ntest], f1{1}'./1e6, c1(:,2:4:end), "GridVisible", "off", "XDisplayLabels", xlabels1,"YDisplayLabels", ylabels, "ColorLimits", [min([c1(:,2:4:end) c2(:,2:4:end)],[],"all"), max([c1(:,2:4:end) c2(:,2:4:end)],[],"all")], "Colormap",jet);
                hs = struct(h);
                xlabel("ntest");
                ylabel("Freq (MHz)");
                ylabel(hs.Colorbar,"dB");
                fig_name="SS Freq Htmp S21 dB at "+string(FFL_or_Phip)+"="+PI_string+" Phi-=" + string(phin_sweep1(end))+ " mA, Right Well "+temp_string+string(dev_string); 
                title({string(fig_name),string(tag)})
                saveas(gcf,string(fig_name) + " " +string(tag)+".fig");
                close(gcf)
                
                figure
                h=heatmap([1:ntest], f2{1}'./1e6, c2(:,2:4:end), "GridVisible", "off", "XDisplayLabels", xlabels2,"YDisplayLabels", ylabels, "ColorLimits", [min([c1(:,2:4:end) c2(:,2:4:end)],[],"all"), max([c1(:,2:4:end) c2(:,2:4:end)],[],"all")], "Colormap",jet);
                hs = struct(h);
                xlabel("ntest");
                ylabel("Freq (MHz)");
                ylabel(hs.Colorbar,"dB");
                fig_name="SS Freq Htmp S21 dB at "+string(FFL_or_Phip)+"="+PI_string+" Phi-=" + string(phin_sweep1(end))+ " mA, Left Well "+temp_string+string(dev_string);   
                title({string(fig_name),string(tag)})
                saveas(gcf,string(fig_name) + " " +string(tag)+".fig");
                close(gcf)
                %Code for heatmaps of S21 Phase vs Phi- & Freq, dont plot
                %{
                figure 
                h=heatmap(flip(phin_sweep1), f1{1}'./1e6, flip(atan(c1(:,4:4:end)./c1(:,3:4:end)),2), "GridVisible", "off", "XDisplayLabels", flip(xlabels1),"YDisplayLabels", ylabels, "ColorLimits", [-1.25,1.25], "Colormap",turbo);
                hs = struct(h);
                xlabel("Phi- (mA)");
                ylabel("Freq (MHz)");
                ylabel(hs.Colorbar,"Phase");
                fig_name="DW Freq Htmp S21 Phase v "+string(FFL_or_Phip)+" & Phi-, Swp1, "+temp_string+string(dev_string)+"loopedx"+string(loop_times);;    
                title({string(fig_name),string(tag)})
                saveas(gcf,string(fig_name) + " " +string(tag)+".fig");
                close(gcf)
            
                figure 
                h=heatmap(phin_sweep2, f2{1}'./1e6, atan(c2(:,4:4:end)./c2(:,3:4:end)), "GridVisible", "off", "XDisplayLabels", xlabels2,"YDisplayLabels", ylabels, "ColorLimits", [-1.25,1.25], "Colormap",turbo);
                hs = struct(h);
                xlabel("Phi- (mA)");
                ylabel("Freq (MHz)");
                ylabel(hs.Colorbar,"Phase");
                fig_name="DW Freq Heatmap S21 Phase v "+string(FFL_or_Phip)+" & Phi-, Swp2, "+temp_string+string(dev_string)+"loopedx"+string(loop_times);;     
                title({string(fig_name),string(tag)})
                saveas(gcf,string(fig_name) + " " +string(tag)+".fig");
                close(gcf)
                %}
            end
        end %plot_flag
    end %test_flag  

   %sweep test for both phip and ffl, use for sweep up and sweep down
   function  [e_test, e_pat_test,  e_ns_test, e_ds_test, e_os_test, Data_test, f_test, dataF_test, dataS_test]= DW_Test(phin_path_arr, phin_sweep_arr)
        if Phin 
            start= phin0;
            sweep=[phin_path_arr]';                                                           %Start location for Phi- with Phi+/FFL in zero state with 1-line Phi-    
            test_sweep=[phin_sweep_arr]';                                           %define Sweep array for actual DW test points along 1-line Phi-
        else
            start= [(phip0 + phin0)./2 , (phin0 - phip0)./2];                       %Start location for Phi- with Phi+ in zero state    
            sweep=[(phip0+phin_path_arr)./2; (phin_path_arr-phip0)./2]';            %define sweep array for initialization Phi- sweep with Phi+ (or FFL) in zero state
            test_sweep=[(phip_PI+phin_sweep_arr)./2; (phin_sweep_arr-phip_PI)./2]'; %define Sweep array for actual DW test points along Phi- with Phi+ (or FFL) in PI state
        end

        for n=1:ntest                                     %For each test point, conduct full initialization sweeps
            n
            c_ma(icb_chan,start);                                               %Set Phip and Phi- to start state
            if FF, c_ma(ffl_chan,ffl0); end  
                                                   %If using FFL, set in 0 state to start
            for i=1:size(sweep,1)                                               %Sweep Phi- with Phi+ or FFL in zero state to Phi- to pivot point
                c_ma(icb_chan,sweep(i,:));     
                %fprintf("Init path Meas. Phi+ = %.4f, Phi- = %.4f\n", c_ma(icb_chan(1))-c_ma(icb_chan(2)), c_ma(icb_chan(1))+c_ma(icb_chan(2)))
            end
            
            if Phin 
                c_ma(icb_chan,phin_path_arr(end));                                                  %1-line Phi- at pivot point
            else
                c_ma(icb_chan,[(phip_PI+phin_path_arr(end))./2; (phin_path_arr(end)-phip_PI)./2]'); %Flip phi+ to PI state with 2-line Phi- at pivot point
            end
            if FF, c_ma(ffl_chan,ffl_PI); end                                                   %If using FFL, set to PI state
            
            for i=1:size(test_sweep,1)                                                                 %Sweep Phi- along Phi+ or FFL = PI up to test point
                c_ma(icb_chan, test_sweep(i,:));
                %fprintf("Test path Meas. Phi+ = %.4f, Phi- = %.4f\n", c_ma(icb_chan(1))-c_ma(icb_chan(2)), c_ma(icb_chan(1))+c_ma(icb_chan(2)))
            end
            
            if Phin
                %fprintf("FFL = %.4f, Phi- = %.4f\n", ffl_PI , phin_sweep_arr(end))
                %fprintf("Test Point Meas. FFL = %.4f, Phi- = %.5f\n", c_ma(ffl_chan), c_ma(icb_chan)) %DO NOT RUN!!! Measuring current kicks JDPD out of state
            else
                if FF
                    %fprintf("FFL = %.2f, Phi- = %.3f\n", ffl_PI , phin_sweep_arr(end))
                    %fprintf("Test Point Meas. FFL = %.4f, Phi- = %.4f\n", c_ma(ffl_chan), c_ma(icb_chan(1))+c_ma(icb_chan(2)))
                else
                    %fprintf("Phi+ = %.2f, Phi- = %.3f\n", phip_PI , phin_sweep_arr(end))
                    %fprintf("Test Point Meas. Phi+ = %.4f, Phi- = %.4f\n", c_ma(icb_chan(1))-c_ma(icb_chan(2)), c_ma(icb_chan(1))+c_ma(icb_chan(2)))
                end
            end

            switch test_flag                                                                 %Run tests with Phi- at test point with Phi+ or FFL = PI in specified order
                case 1
                    [e_test(n), e_pat_test(n), e_ns_test(n), e_ds_test(n), e_os_test(n)]=test_ER();  
                case 2
                    [Data_test{n}, f_test{n}]=test_VNA();
                case 3
                    [Data_test{n}, f_test{n}]=test_VNA();
                    [e_test(n), e_pat_test(n), e_ns_test(n), e_ds_test(n), e_os_test(n)]=test_ER(); 
                case 4
                    [e_test(n), e_pat_test(n), e_ns_test(n), e_ds_test(n), e_os_test(n)]=test_ER();  
                    [Data_test{n}, f_test{n}]=test_VNA();
           end
        end
        if test_flag==1
            Data_test= {nan};
            f_test= {nan};
        elseif test_flag==2
            e_test=nan;
            e_pat_test=nan;   
            e_ns_test=nan;  
            e_ds_test=nan;  
            e_os_test=nan;  
        end
    end

    function [e_bit_arg, e_pat_arg, e_ns_arg, e_ds_arg, e_os_arg]=test_ER()     
        %test pattern single shot, one pattern
        [p, b_er, ns_er, ds_er, os_er] = pattern_test_spliterrors(pattest,edges,[],0);  %run test pattern, no patclears, assume trigger on variable "edges" (defined as 1 or rising edge), split and categorize errors
           
        e_pat_arg= 1-p;
        e_bit_arg= b_er/(bits);
        e_ns_arg=  ns_er/(bits);
        e_ds_arg=  ds_er/(bits);
        e_os_arg=  os_er/(bits);
        %%%fprintf("Pat ER is %.2f and Bit ER rate is %.4f [%.2f NS %.2f DS %.2f OS]  \n", e_pat_arg, e_bit_arg, e_ns_arg, e_ds_arg, e_os_arg);
    end  

    function [Data_arg, f_arg]=test_VNA()
        [f_arg, dataF_arg, dataS_arg]= v.getData;
        Data_arg=[f_arg' dataF_arg' dataS_arg(1:2:end)'  dataS_arg(2:2:end)']; 
        avgp=mean(Data_arg(idx1:idx2,2));
        %%%fprintf("@%.2f GHz, Avg S21 = %.2f dB\n", round(mean(sample_range(1,:))./1e9,2), avgp)       
    end
end