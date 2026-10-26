%*********************************************Double Well VNA ER Updated***************************************************

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
    b1=[];
    b2=[];
    e1_pat=[];
    e2_pat=[];

    e1_ns=[]; 
    e1_ds=[];
    e1_os=[];
    e2_ns=[];
    e2_ds=[];
    e2_os=[];

    e1_iter=[];
    e1_pat_iter =[];
    e1_ns_iter =[];
    e1_ds_iter =[];
    e1_os_iter =[]; 

    e2_iter=[];
    e2_pat_iter =[];
    e2_ns_iter =[];
    e2_ds_iter =[];
    e2_os_iter =[];
    
    edges;
    splitcurrent = 0;
    if test_flag~=2                                             %If bitwise is true, assign bit value based on test pattern expected 1's
        if islogical(bitwise) && bitwise
            bit_flag=true;
            bitwise=sum(PATTERN(pattest).vdat,"all");
        elseif (islogical(bitwise) && ~bitwise) || bitwise<=0   %If bitwise is false or 0, still assign bit value for full bitwise tests but later only plot pattern error rate
            bit_flag=false;
            bitwise=sum(PATTERN(pattest).vdat,"all");
        elseif bitwise>0                                        %If bitwise is a positive whole number, use THAT number for BER calculations
            bit_flag=true;
            bitwise=round(bitwise);
        end
     end

    %VNA initializing
    f1={};
    dataF1={};
    dataS1={};
    Data1={};
    AvgPwr1={};
    f2={};
    dataF2={};
    dataS2={};
    Data2={};
    AvgPwr2={};   
    avgp=[];
    Data1_iter={}; 
    f1_iter ={}; 
    dataF1_iter={}; 
    dataS1_iter={}; 
    Data2_iter={}; 
    f2_iter ={}; 
    dataF2_iter={}; 
    dataS2_iter={}; 
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
        if isempty(idx1) || isnan(idx1), idx1=1; end
        if isempty(idx2) || isnan(idx2), idx2=length(rang); end
    end

    % icb is array of [DCFlux1, DCFlux2], also called icb(1) and icb(2), they are currents applied to induce flux in JJ1 and JJ2. 
    % Phi+= DCFlux1 - DCFlux2, Phi-= DCFlux1 + DCFlux2
    % DCFlux1 = (Phi- + Phi+)/2, DCFlux2 = (Phi- - Phi+)/2

    phin0=phin_path1(1);                    %Starting Phi minus for both sweeps
    phin_path2=2*phin_path1(1)-phin_path1;  %Define return initialization sweep based on first initialization sweep
    phin_sweep2=2*phin_sweep1(end)-phin_sweep1;%Define return test sweep based on first test sweep

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

    % if loop_times==1            %Single iteration of full DW test
        disp("************************************************************Sweep1")
        [e1, e1_pat, e1_ns, e1_ds, e1_os, Data1,  f1, dataF1, dataS1]=DW_Test(phin_path1,phin_sweep1);    
        disp("************************************Sweep2")
        [e2, e2_pat, e2_ns, e2_ds, e2_os, Data2,  f2, dataF2, dataS2]=DW_Test(phin_path2,phin_sweep2); 
    %     e1=mean(e1);
    %     e2=mean(e2);
    %     e1_pat=mean(e1_pat);
    %     e2_pat=mean(e2_pat);
    %     e1_ns=mean(e1_ns);
    %     e2_ns=mean(e2_ns);
    %     e1_ds=mean(e1_ds);
    %     e2_ds=mean(e2_ds);
    %     e1_os=mean(e1_os);
    %     e2_os=mean(e2_os);    

    %elseif loop_times>1         %Do full DW tests iterated loop_times, average all results
        
            % fprintf("*****************************Iteration %d Sweep1*******************\n",loops)
            % [e1_iter(:),e1_pat_iter(:), e1_ns_iter(:), e1_ds_iter(:), e1_os_iter(:),  Data1_iter(:), f1_iter(:), dataF1_iter(:), dataS1_iter(:)]=DW_Test(phin_path1,phin_sweep1);         
            % fprintf("*****************Iteration %d Sweep2*****\n",loops)
            % [e2_iter(loops,:),e2_pat_iter(loops,:), e2_ns_iter(loops,:), e2_ds_iter(loops,:), e2_os_iter(loops,:),  Data2_iter(loops,:), f2_iter(loops,:), dataF2_iter(loops,:), dataS2_iter(loops,:)][...]
            % 
        %Average data for each loop independently (combine all results from each ntest)
        % e1=mean(e1_iter,2);
        % e2=mean(e2_iter,2);
        % e1_pat=mean(e1_pat_iter,2);
        % e2_pat=mean(e2_pat_iter,2);
        % e1_ns=mean(e1_ns_iter,2);
        % e2_ns=mean(e2_ns_iter,2);
        % e1_ds=mean(e1_ds_iter,2);
        % e2_ds=mean(e2_ds_iter,2);
        % e1_os=mean(e1_os_iter,2);
        % e2_os=mean(e2_os_iter,2);

        % e1=e1_iter;
        % e2=e2_iter;
        % e1_pat=e1_pat_iter;
        % e2_pat=e2_pat_iter;
        % e1_ns=e1_ns_iter;
        % e2_ns=e2_ns_iter;
        % e1_ds=e1_ds_iter;
        % e2_ds=e2_ds_iter;
        % e1_os=e1_os_iter;
        % e2_os=e2_os_iter;




%{
        Data1 = cell(1, size(Data1_iter, 2));
        f1 = cell(1, size(f1_iter, 2));
        dataF1 = cell(1, size(dataF1_iter, 2));
        dataS1 = cell(1, size(dataS1_iter, 2));
        Data2 = cell(1, size(Data2_iter, 2));
        f2 = cell(1, size(f2_iter, 2));
        dataF2 = cell(1, size(dataF2_iter, 2));
        dataS2 = cell(1, size(dataS2_iter, 2));        
        for col = 1:size(Data1_iter, 2)
            Data1{col} = mean(cat(3, Data1_iter{:, col}), 3, 'omitnan');
            f1{col} = mean(cat(3, f1_iter{:, col}), 3, 'omitnan');
            dataF1{col} = mean(cat(3, dataF1_iter{:, col}), 3, 'omitnan');
            dataS1{col} = mean(cat(3, dataS1_iter{:, col}), 3, 'omitnan');        
            Data2{col} = mean(cat(3, Data2_iter{:, col}), 3, 'omitnan');
            f2{col} = mean(cat(3, f2_iter{:, col}), 3, 'omitnan');
            dataF2{col} = mean(cat(3, dataF2_iter{:, col}), 3, 'omitnan');
            dataS2{col} = mean(cat(3, dataS2_iter{:, col}), 3, 'omitnan');
        end
    end
%}

    c_ma(icb_chan,[0,0]);                                               %turn off Phi+, Phi-, and FFL
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
    if ~isempty(und_idx)
        dev_string=string(icb_char(1:und_idx-1));
    else
        dev_string="";
    end
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
            % Code for plot of DW comp Total bit ER and also ER due to double switch, offcycle switch               
            if bit_flag
                figure
                hold on
                leg1="Right Well "+string(mean(e1)*100)+"% Tot BER";
                leg2="Left Well "+ string(mean(e2)*100)+"% Tot BER";
                leg3="Right Well "+string(mean(e1_ds+e1_os)*100)+"% DS+OS BER";
                leg4="Left Well "+string(mean(e2_ds+e2_os)*100)+"% DS+OS BER";
    
                plot([1:ntest],e1,"o","DisplayName", string(leg1));
                plot([1:ntest],e2,'*',"DisplayName", string(leg2));
                plot([1:ntest],e1_ds+e1_os,"DisplayName", string(leg3), "Visible", "off", "YLimInclude","off","HandleVisibility","off") 
                plot([1:ntest],e2_ds+e2_os,"DisplayName", string(leg4), "Visible", "off", "YLimInclude","off","HandleVisibility","off")  
                ylim([0,1]);
                xlabel("ntest")
                ylabel("BER")
                legend
                fig_name="SS Bit ER at "+string(FFL_or_Phip)+"="+PI_string+" Phi-=" + string(phin_sweep1(end))+ " mA, "+ temp_string+string(dev_string)+string(ntest)+"x"+string(bitwise)+"b";    
                title({string(fig_name),string(tag)})
                saveas(gcf,string(fig_name) + " " +string(tag)+".fig");
                hold off
                close(gcf)
                
                figure
                histogram(e1,50,"FaceColor","blue","BinLimits", [0 1])
                xlabel("BER by Each Test Pattern")
                ylabel("Count")
                fig_name="Hist SS Right Well Pat BER, "+string(FFL_or_Phip)+"="+PI_string+" Phi-=" + string(phin_sweep1(end))+ " mA, "+ temp_string+string(dev_string)+string(ntest)+"x"+string(bitwise)+"b";
                title({string(fig_name),string(tag)})
                xlim([0,1.05])
                saveas(gcf,string(fig_name) + " " +string(tag)+".fig");
                close(gcf)

                figure
                histogram(e2,50,"FaceColor","red","BinLimits", [0 1])
                xlabel("BER by Each Test Pattern")
                ylabel("Count")
                fig_name="Hist SS Left Well Pat BER "+string(FFL_or_Phip)+"="+PI_string+" Phi-=" + string(phin_sweep1(end))+ " mA, "+ temp_string+string(dev_string)+string(ntest)+"x"+string(bitwise)+"b";
                title({string(fig_name),string(tag)})
                xlim([0,1.05])
                saveas(gcf,string(fig_name) + " " +string(tag)+".fig");
                close(gcf)

            elseif ~bit_flag
              % Code for plot of DW comp Pattern ER
                figure
                hold on
                leg1="Right Well "+string(mean(e1_pat)*100)+"% Pat ER";
                leg2="Left Well "+string(mean(e2_pat)*100)+"% Pat ER";
    
                plot([1:ntest],e1_pat,"o","DisplayName", string(leg1));
                plot([1:ntest],e2_pat,'*',"DisplayName", string(leg2));
                ylim([0,1]);
                xlabel("ntest")
                ylabel("Pat ER")
                legend
                fig_name="SS Pat ER at "+string(FFL_or_Phip)+"="+PI_string+" Phi-=" + string(phin_sweep1(end))+ " mA, "+ temp_string+string(dev_string)+string(ntest)+"x"+string(bitwise)+"b";     
                title({string(fig_name),string(tag)})
                saveas(gcf,string(fig_name) + " " +string(tag)+".fig");
                hold off
                close(gcf)
            end %bitflag
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
                plot([1:ntest], p1(:,2:4:end),"o","DisplayName", string(leg3)+string(round(mean(sample_range(k,:))./1e9,2)) + "GHz", "Visible", "off", "YLimInclude","off","HandleVisibility","off")    [...]
                plot([1:ntest], p2(:,2:4:end),'*',"DisplayName", string(leg4)+string(round(mean(sample_range(k,:))./1e9,2)) + "GHz", "Visible", "off", "YLimInclude","off","HandleVisibility","off") 
                plot([1:ntest], p1(:,3:4:end),"o","DisplayName", string(leg5)+string(round(mean(sample_range(k,:))./1e9,2)) + "GHz", "Visible", "off", "YLimInclude","off","HandleVisibility","off") 
                plot([1:ntest], p2(:,3:4:end),'*',"DisplayName", string(leg6)+string(round(mean(sample_range(k,:))./1e9,2)) + "GHz", "Visible", "off", "YLimInclude","off","HandleVisibility","off") 
                plot([1:ntest], p1(:,4:4:end),"o","DisplayName", string(leg7)+string(round(mean(sample_range(k,:))./1e9,2)) + "GHz", "Visible", "off", "YLimInclude","off","HandleVisibility","off") 
                plot([1:ntest], p2(:,4:4:end),'*',"DisplayName", string(leg8)+string(round(mean(sample_range(k,:))./1e9,2)) + "GHz", "Visible", "off", "YLimInclude","off","HandleVisibility","off")    [...]
            end
   
            axis tight
            hold off
            xlabel("ntest")
            ylabel("S21")
            legend
            fig_name="SS S21 dB at "+string(FFL_or_Phip)+"="+PI_string+" Phi-=" + string(phin_sweep1(end))+ " mA, "+temp_string+string(dev_string)+string(round(mean(sample_range(1,:))./1e9,2)) + "GHz"[...]
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
                h=heatmap([1:ntest], f1{1}'./1e6, c1(:,2:4:end), "GridVisible", "off", "XDisplayLabels", xlabels1,"YDisplayLabels", ylabels, "ColorLimits", [min([c1(:,2:4:end) c2(:,2:4:end)],[],"all") [...]
                hs = struct(h);
                xlabel("ntest");
                ylabel("Freq (MHz)");
                ylabel(hs.Colorbar,"dB");
                fig_name="SS Freq Htmp S21 dB at "+string(FFL_or_Phip)+"="+PI_string+" Phi-=" + string(phin_sweep1(end))+ " mA, Right Well "+temp_string+string(dev_string); 
                title({string(fig_name),string(tag)})
                saveas(gcf,string(fig_name) + " " +string(tag)+".fig");
                close(gcf)
                
                figure
                h=heatmap([1:ntest], f2{1}'./1e6, c2(:,2:4:end), "GridVisible", "off", "XDisplayLabels", xlabels2,"YDisplayLabels", ylabels, "ColorLimits", [min([c1(:,2:4:end) c2(:,2:4:end)],[],"all") [...]
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
                h=heatmap(flip(phin_sweep1), f1{1}'./1e6, flip(atan(c1(:,4:4:end)./c1(:,3:4:end)),2), "GridVisible", "off", "XDisplayLabels", flip(xlabels1),"YDisplayLabels", ylabels, "ColorLimits", [[...]
                hs = struct(h);
                xlabel("Phi- (mA)");
                ylabel("Freq (MHz)");
                ylabel(hs.Colorbar,"Phase");
                fig_name="DW Freq Htmp S21 Phase v "+string(FFL_or_Phip)+" & Phi-, Swp1, "+temp_string+string(dev_string)+"loopedx"+string(loop_times);;    
                title({string(fig_name),string(tag)})
                saveas(gcf,string(fig_name) + " " +string(tag)+".fig");
                close(gcf)
            
                figure 
                h=heatmap(phin_sweep2, f2{1}'./1e6, atan(c2(:,4:4:end)./c2(:,3:4:end)), "GridVisible", "off", "XDisplayLabels", xlabels2,"YDisplayLabels", ylabels, "ColorLimits", [-1.25,1.25], "Colorm[...]
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
        e_test=[];      %Temporary function variable initialization 
        e_pat_test=[];
        e_ns_test=[];
        e_ds_test=[];
        e_os_test=[];
        f_test={};
        dataF_test={};
        dataS_test={};
        Data_test={};
        
        start= [(phip0 + phin0)./2 , (phin0 - phip0)./2];                       %Start location for Phi- with Phi+ in zero state    
        sweep=[(phip0+phin_path_arr)./2; (phin_path_arr-phip0)./2]';            %define sweep array for initialization Phi- sweep with Phi+ (or FFL) in zero state
        test_sweep=[(phip_PI+phin_sweep_arr)./2; (phin_sweep_arr-phip_PI)./2]'; %define Sweep array for actual DW test points along Phi- with Phi+ (or FFL) in PI state
      
        for n=1:ntest                                     %For each test point, conduct full initialization sweeps
            c_ma(icb_chan,start);                                               %Set Phip and Phi- to start state
            if FF, c_ma(ffl_chan,ffl0); end                                     %If using FFL, set in 0 state to start
            for i=1:size(sweep,1)                                               %Sweep Phi- with Phi+ or FFL in zero state to Phi- to pivot point
                c_ma(icb_chan,sweep(i,:));     
                %fprintf("Init path Meas. Phi+ = %.4f, Phi- = %.4f\n", c_ma(icb_chan(1))-c_ma(icb_chan(2)), c_ma(icb_chan(1))+c_ma(icb_chan(2)))
            end
            
            c_ma(icb_chan,[(phip_PI+phin_path_arr(end))./2; (phin_path_arr(end)-phip_PI)./2]');   %Flip phi+ to PI state with phi- at pivot point
            if FF, c_ma(ffl_chan,ffl_PI); end                                                   %If using FFL, set to PI state
            
            for i=1:size(test_sweep,1)                                                                 %Sweep Phi- along Phi+ or FFL = PI up to test point
                c_ma(icb_chan, test_sweep(i,:));
                %fprintf("Test path Meas. Phi+ = %.4f, Phi- = %.4f\n", c_ma(icb_chan(1))-c_ma(icb_chan(2)), c_ma(icb_chan(1))+c_ma(icb_chan(2)))
            end           
            if FF
                %fprintf("FFL = %.2f, Phi- = %.3f\n", ffl_PI , phin_sweep_arr(end))
                %fprintf("Test Point Meas. FFL = %.4f, Phi- = %.4f\n", c_ma(ffl_chan), c_ma(icb_chan(1))+c_ma(icb_chan(2)))
            else
                %fprintf("Phi+ = %.2f, Phi- = %.3f\n", phip_PI , phin_sweep_arr(end))
                %fprintf("Test Point Meas. Phi+ = %.4f, Phi- = %.4f\n", c_ma(icb_chan(1))-c_ma(icb_chan(2)), c_ma(icb_chan(1))+c_ma(icb_chan(2)))
            end

            switch test_flag                                                                 %Run tests with Phi- at test point with Phi+ or FFL = PI in specified order
                case 1
                    [e_test(n), e_pat_test(n), e_ns_test(n), e_ds_test(n), e_os_test(n)]=test_ER();  
                case 2
                    [Data_test{n}, f_test{n}, dataF_test{n}, dataS_test{n}]=test_VNA();
                case 3
                    [Data_test{n}, f_test{n}, dataF_test{n}, dataS_test{n}]=test_VNA();
                    [e_test(n), e_pat_test(n), e_ns_test(n), e_ds_test(n), e_os_test(n)]=test_ER(); 
                case 4
                    [e_test(n), e_pat_test(n), e_ns_test(n), e_ds_test(n), e_os_test(n)]=test_ER();  
                    [Data_test{n}, f_test{n}, dataF_test{n}, dataS_test{n}]=test_VNA();
           end
        end
        if test_flag==1
            Data_test= {nan};
            f_test= {nan};
            dataF_test= {nan};
            dataS_test= {nan};
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
        [p, b_er, ns_er, ds_er, os_er] = pattern_test_spliterrors(pattest,edges,[],0);  %run test pattern, no patclears, assume trigger on variable "edges" (defined as 1 or rising edge), split and cat[...]
           
        e_pat_arg= 1-p;
        e_bit_arg= b_er/(bitwise);
        e_ns_arg=  ns_er/(bitwise);
        e_ds_arg=  ds_er/(bitwise);
        e_os_arg=  os_er/(bitwise);
        %%%fprintf("Pat ER is %.2f and Bit ER rate is %.4f [%.2f NS %.2f DS %.2f OS]  \n", e_pat_arg, e_bit_arg, e_ns_arg, e_ds_arg, e_os_arg);
    end  

    function [Data_arg, f_arg, dataF_arg, dataS_arg]=test_VNA()
        [f_arg, dataF_arg, dataS_arg]= v.getData;
        Data_arg=[f_arg' dataF_arg' dataS_arg(1:2:end)'  dataS_arg(2:2:end)']; 
        avgp=mean(Data_arg(idx1:idx2,2));
        %%%fprintf("@%.2f GHz, Avg S21 = %.2f dB\n", round(mean(sample_range(1,:))./1e9,2), avgp)       
    end


%Old version of code

%{
    function  [e_test, e_pat_test, Data_test, f_test, dataF_test, dataS_test]= Sweep_Phip(phin_path_arr, phin_sweep_arr)
        e_test=[];
        e_pat_test=[];
        f_test={};
        dataF_test={};
        dataS_test={};
        Data_test={};

        start= [(phip0 + phin0)./2 , (phin0 - phip0)./2]; %Start location for Phi- with Phi+=0
        icb=start;
        test_sweep=[(phip_PI+phin_sweep_arr)./2; (phin_sweep_arr-phip_PI)./2]'; %Sweep array for first actual DW test points along Phi-
        for j=1:length(test_sweep)                                              %For each test point, conduct full initialization sweeps
            c_ma(icb_chan,start);
            sweep=[(phip0+phin_path_arr)./2; (phin_path_arr-phip0)./2]';        %sweep array for first initial Phi- sweep along Phi+ = 0
            for i=1:length(sweep)                                               %Sweep Phi- along Phi+=0 to pivot point
                c_ma(icb_chan,sweep(i,:));           
            end
            %c_ma(icb_chan,[(phip0+phin_path_var(end))./2; (phin_path_var(end)-phip0)./2]');                                %phi+ in 0 state with phi- at pivot point
            c_ma(icb_chan,[(phip_PI+phin_path_arr(end))./2; (phin_path_arr(end)-phip_PI)./2]');   %Flip phi+ to PI state with phi- at pivot point
            for i=1:j                                                                           %Sweep Phi- along Phi+ = PI up to test point, first sweep
                c_ma(icb_chan, test_sweep(i,:));
            end           
            icb= test_sweep(j,:);                                                            %Run tests with Phi- at test point with Phi+=PI 
            fprintf("Phi+ = %.2f, Phi- = %.3f\n", phip_PI , phin_sweep_arr(j))            
            switch test_flag
                case 1
                    [e_test(j), e_pat_test(j)]=test_ER();  
                case 2
                    [Data_test{j}, f_test{j}, dataF_test{j}, dataS_test{j}]=test_VNA();
                case 3
                    [Data_test{j}, f_test{j}, dataF_test{j}, dataS_test{j}]=test_VNA();
                    [e_test(j), e_pat_test(j)]=test_ER();  
                case 4
                    [e_test(j), e_pat_test(j)]=test_ER();  
                    [Data_test{j}, f_test{j}, dataF_test{j}, dataS_test{j}]=test_VNA();
           end
        end
    end
  
    function [e_test, e_pat_test, Data_test, f_test, dataF_test, dataS_test]=Sweep_FFL(phin_path_arr, phin_sweep_arr)
        e_test=[];
        e_pat_test=[];
        f_test={};
        dataF_test={};
        dataS_test={};
        Data_test={};

        start= [phin0./2 , phin0./2];
        icb=start;
        test_sweep=[phin_sweep_arr./2; phin_sweep_arr./2]';
        c_ma(icb_chan,start);
        for j=1:length(test_sweep)        
            c_ma(ffl_chan,phip0);
            sweep=[phin_path_arr./2; phin_path_arr./2]';
            for i=1:length(sweep)
                c_ma(icb_chan,sweep(i,:));
            end
            %c_ma(ffl_chan,phip0); %phi plus in 0 state with phi minus at end of inital path
            c_ma(ffl_chan,phip_PI); % phi plus in PI state with phi minus at end of inital path   
            for i=1:j
                c_ma(icb_chan, test_sweep(i,:));
            end
            icb= test_sweep(j,:);
            fprintf("FFL = %.2f, Phi- = %.3f\n", phip_PI , phin_sweep_arr(j))
            switch test_flag
                case 1
                    [e_test(j), e_pat_test(j)]=test_ER();  
                case 2
                    [Data_test{j}, f_test{j}, dataF_test{j}, dataS_test{j}]=test_VNA();
                case 3
                    [Data_test{j}, f_test{j}, dataF_test{j}, dataS_test{j}]=test_VNA();
                    [e_test(j), e_pat_test(j)]=test_ER();  
                case 4
                    [e_test(j), e_pat_test(j)]=test_ER();  
                    [Data_test{j}, f_test{j}, dataF_test{j}, dataS_test{j}]=test_VNA();
            end
        end
    end
%}

end
