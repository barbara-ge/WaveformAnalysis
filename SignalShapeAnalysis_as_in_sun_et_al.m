%%calculatine positive spoikes vs negative spikes

clear all
close all
folderList = [{'week2'},{'div19'},{'week4'}];
populationList = [{'NS'},{'8020'}, {'5050'}];
cd('/Volumes/Elements/spring2022/final/')
load('idx_cell.mat')
SpikeFeature = {};
SpikeShapePop = {};
TotCountSpike = {};
FeaturePopulation = {};
AmplitudePopulation = {};

for ff = 1: length(folderList)
         ff
    %cd(['' folderList{1,ff} '/CSV/'])
    cd(['' folderList{1,ff} '/'])
    
%     if ff == 2
%         start = 3;
%         endList = 3;
%     else
%         start = 1;
%         endList = 2;
%     end

      
    for pp = 1: length(populationList)
        pp

%         if pp == 3

            cd(['' populationList{1,pp} ''])
%         else
%             cd(['' populationList{1,pp} '/analysed/'])
%         end

        mea_folder = dir('*');
        posMEAPop= [];
        tsMEAPop = [];
        csMEAPop = [];
        fsMEAPop = [];
        rsMEAPop = [];
        negMEAPop = [];
        otherMEAPop = [];
        slopeMEAPop = [];
        durationMEAPop = [];
        Amp1Pop = [];
        Amp2Pop = [];
        AmpPosPop = [];
        AmpNegPop = [];

        for ii= 3:length(mea_folder)
            mea_act_mat = zeros(length(mea_folder),2);

            mea_num = mea_folder(ii).name;
            mea_num
            spikecounts = {};
            SpikeShape = {};
            posMEA = [];
            tsMEA = [];
            csMEA = [];
            fsMEA = [];
            rsMEA = [];
            negMEA = [];
            otherMEA = [];
            slopeMEA = [];
            durationMEA = [];
            AmpPosMEA = [];
            AmpNegMEA = [];
            amp1MEA = [];
            amp2MEA = [];

            if isdir(mea_num)==0 || isempty(mea_num) == 1 
                continue
            else

                path_cd = ['' mea_num '/'];
                
                cd(path_cd)

                SpikeTypeMat = [];
                SlopeMat = [];
                DurationMat = [];
                ampPosMat = [];
                ampNegMat = [];
                Amp1Mat = [];
                Amp2Mat = [];

                for cc = 1 : 60
                    pos= 0;
                    rs = 0;
                    fs = 0;
                    ts = 0;
                    neg = 0;
                    cs = 0;
                    other = 0;

                    numCC = IDs{cc,1};

                    if strcmp(numCC,'Ref') == 1
                        continue
                    else
                        filenameCh = ['times_ch' numCC 'data.mat'];

                        if exist(filenameCh)

                            load(filenameCh)
                            clusterNum = unique(cluster_class(:,1));

                            for cn = 1: length(clusterNum)

                                aa = clusterNum(cn);
                                signalClustIdx = find(cluster_class(:,1)==aa);

                                if length(signalClustIdx)==1
                                    meanSignal = spikes(signalClustIdx,:);
                                else
                                    spikes_cluster = spikes(signalClustIdx,:);

                                    meanSignal = mean(spikes_cluster(:,:));
                                end

                                spikeAmp = meanSignal(1,20);
                                baseline = mean([meanSignal(1,1),meanSignal(1,end)]);

                                if  spikeAmp > 0
                                    spikeType = 1;
                                    pos = pos+length(signalClustIdx);
                                    [pks, index] = findpeaks(-meanSignal);
                                    a = find(index>20);
                                    b = find(index < 20);
                                    ampPosMat = [ampPosMat; spikeAmp];

                                    if  isempty(b) == 1
                                        IndexPreviousPeakCorrected = 1;
                                        amp1 = meanSignal(1,20) - meanSignal(IndexPreviousPeakCorrected);
                                        slope1 = amp1 ./ 20 ;
                                    else
                                        IndexPreviousPeakCorrected = index(b(end));
                                        amp1 = spikeAmp - baseline; %meanSignal(1,1); %spike - baseline
                                        
%                                         if -pks(b(end)) > 0
%                                             slope1 = meanSignal(1,20) ./ 20 ;
%                                         else
                                        slope1 = (meanSignal(1,20) - (-pks(b(end)))) ./ (20 - index(b(end)));
%                                         end
                                
                                    end

                                    if  isempty(a) == 1
                                        [minValue,closestIndex] = min(abs(meanSignal(20:end)-meanSignal(1,1)));
                                        amp2 =  minValue - baseline; %- meanSignal(1,1);
                                        Duration = closestIndex*0.04;
                                        %PeakToPeak = ((closestIndex+20) - IndexPreviousPeakCorrected)*0.04;
                                    else
                                        amp2 = (-pks(a(1))) - baseline; %- meanSignal(1,1);
                                       % PeakToPeak = (index(a(1)) - IndexPreviousPeakCorrected)*0.04;

                                        Duration = (index(a(1)) - 20)*0.04; %from time points to ms
                                        [minValue,closestIndex] = min(abs(meanSignal(index(a(1)):end)-meanSignal(1,1)));
                                        slope2 = (minValue - amp2) ./ (closestIndex - index(a(1)));

                                    end
                                else
                                    ampNegMat = [ampNegMat; spikeAmp];
                                    
                                    neg = neg + length(signalClustIdx);
                                    [pks, index] = findpeaks(meanSignal);
                                    %  Signal2 = meanSignal(1,index(1):index(2));

                                    %[maximum2 lmax2] = max(meanSignal(1,index>lmin));
                                    a = find(index>20);
                                    b = find(index < 20);

                                    if  isempty(b) == 1
                                        IndexPreviousPeakCorrected = 1;
                                        amp1 = meanSignal(1,20) - meanSignal(IndexPreviousPeakCorrected);
                                        slope1 = amp1 ./ 20 ;
                                        FirstPeakThroughR = 0;
                                    else
                                        IndexPreviousPeakCorrected = index(b(end));
                                        amp1 = spikeAmp - baseline; %meanSignal(1,1); %spike - baseline
                                        FirstPeak = pks(b(end)) - baseline; %meanSignal(1,1); %1st peak - baseline
                                        slope1 = (meanSignal(1,20) - pks(b(end))) ./ (20 - index(b(end)));
                                        FirstPeakThroughR = abs(FirstPeak/amp1);

                                    end

                                    if  isempty(a) == 1
                                        [minValue,closestIndex] = min(abs(meanSignal(20:end)- baseline)); %-meanSignal(1,1)));
                                        IndexLastPeakCorrected = 60;
                                        amp2 = minValue - baseline; %- meanSignal(1,1);
                                        Duration = closestIndex*0.04;
                                        %PeakToPeak = ((closestIndex+20) - index(b(end)))*0.04;
                                    else
                                        amp2 = pks(a(1)) - baseline; %- meanSignal(1,1);
                                        IndexLastPeakCorrected = index(a(1));
                                        %PeakToPeak = (index(a(1)) - index(b(end)))*0.04;
                                        Duration = (index(a(1)) - 20)*0.04; %from time points to ms
                                        [minValue,closestIndex] = min(abs(meanSignal(index(a(1)):end)- baseline)); %-meanSignal(1,1)));
                                        slope2 = (minValue - amp2) ./ (closestIndex - index(a(1)));

                                    end


                                    if FirstPeakThroughR > 0.1
                                        
                                        PeakToPeak = (IndexLastPeakCorrected - IndexPreviousPeakCorrected)*0.04;
                                        if PeakToPeak > 1
                                            cs = cs +length(signalClustIdx);
                                            spikeType = 2;
                                            other = other + length(signalClustIdx);
                                        else
                                            ts = ts +length(signalClustIdx);
                                            other = other + length(signalClustIdx);
                                            spikeType = 3;
                                        end
                                        %[minimum2 lmin2] = min(meanSignal(1,index(a(1)) : end));
                                        %indexMinCorrected = index(a(1))+lmin2;

                                    else
                                        if Duration < 0.3
                                            fs = fs +length(signalClustIdx);
                                            spikeType = 4;
                                        else
                                            rs = rs +length(signalClustIdx);
                                            spikeType = 5;
                                        end

                                    end
                                end

                                %                         spikePoints(cc, 1) = amp1;
                                %                         spikePoints(cc, 2) = amp2;
                                %                         spikePoints(cc, 3) = SpikeDuration;
                                %
                                %
                                %
                                %
                                %                         spikePoints(cc, 4) = slope1;
                                %
                                %                         spikePoints(cc,5) = slope2;
                                SpikeTypeMat = [SpikeTypeMat; spikeType];
                                SlopeMat = [SlopeMat; slope1];
                                DurationMat = [DurationMat; Duration];
                                Amp1Mat = [Amp1Mat; amp1];
                                Amp2Mat = [Amp2Mat; amp2];
                            end

                            spikecounts{cc}(1,1) = pos;
                            spikecounts{cc}(1,2) = cs;
                            spikecounts{cc}(1,3) = ts;
                            spikecounts{cc}(1,4) = fs;
                            spikecounts{cc}(1,5) = rs;
                            spikecounts{cc}(1,6) = neg;
                            spikecounts{cc}(1,7) = other;
                            SpikeShape{cc} = SpikeTypeMat;
                             

                        else

                            continue

                        end

                    end
                    posMEA = [posMEA; pos/length(spikes(:,:))];
                    tsMEA = [tsMEA; ts/length(spikes(:,:))];
                    csMEA = [csMEA; cs/length(spikes(:,:))];
                    fsMEA = [fsMEA; fs/length(spikes(:,:))];
                    rsMEA = [rsMEA; rs/length(spikes(:,:))];
                    negMEA = [negMEA; neg/length(spikes(:,:))];
                    otherMEA = [otherMEA; other/length(spikes(:,:))];
                    AmpPosMEA = [AmpPosMEA; ampPosMat];
                    AmpNegMEA = [AmpNegMEA; ampNegMat];
                    amp1MEA = [amp1MEA; Amp1Mat];
                    amp2MEA = [amp2MEA; Amp2Mat];
                    slopeMEA = [slopeMEA; SlopeMat];
                    durationMEA = [durationMEA; DurationMat];
                end

                posMEAPop = [posMEAPop; posMEA];
                tsMEAPop = [tsMEAPop; tsMEA];
                csMEAPop = [csMEAPop; csMEA];
                fsMEAPop = [fsMEAPop; fsMEA];
                rsMEAPop = [rsMEAPop; rsMEA];
                negMEAPop = [negMEAPop; negMEA];
                otherMEAPop = [otherMEAPop; otherMEA];
                slopeMEAPop = [slopeMEAPop; slopeMEA];
                durationMEAPop = [durationMEAPop; durationMEA];
                AmpPosPop = [AmpPosPop; AmpPosMEA];
                AmpNegPop = [AmpNegPop; AmpNegMEA];
                Amp1Pop = [Amp1Pop; amp1MEA];
                Amp2Pop = [Amp2Pop; amp2MEA];
                SpikeFeature{ii,pp} = spikecounts;
                SpikeShapePop{ii,pp} = SpikeShape;


                cd ..
            end
        end
        
        TotCountSpikes{1,pp} = {posMEAPop, tsMEAPop, csMEAPop, fsMEAPop, rsMEAPop, negMEAPop, otherMEAPop};
        FeaturePopulation{1,pp} = {slopeMEAPop, durationMEAPop};
        AmplitudePopulation{1,pp} = {AmpPosPop, AmpNegPop, Amp1Pop, Amp2Pop};
    
        cd ..
    end
    cd ..
    save(['SpikeShapeAnalysis_final_peak_to_peak_' folderList{1,ff} '.mat'], 'SpikeFeature', 'SpikeShapePop', 'TotCountSpikes', 'FeaturePopulation', 'AmplitudePopulation')

    
end





