function [ res ] = makeViralProbesAlternatingAdapters( sequenceFileName, outFileName, ...
        adapterA_5, adapterA_3, adapterB_5, adapterB_3, overlap, fragLength, rc)
    
    %{
        Produces a FASTA file of probes with alternating A and B adapters
        on the 5' and 3' side, as well as the reverse complement if specified. 
        Portions of the genomes with Ns are excluded but a probe is added
        to bump right up against the Ns
        
        Takes the following inputs:
        sequenceFileName: FASTA file containing the genomes being tiled
        outFileName: Where to write the file
        adapterA_5: 5' sequence of the A adapter
        adapterA_3: 3' sequence of the A adapter
        adapterB_5: 5' sequence of the B adapter
        adapterB_3: 3' sequence of the B adapter
        overlap: the number of bases that the probes should overlap 
        fragLength: the length of the probe 
        rc: set to 1 if you want reverse complement probes
        
        Example:
                
        rc=1;
        overlap=50;
        fragLength=100;

        sequenceFileName=strcat('C:\Temp\ViralGenome.fasta');
        outFileName=strcat('C:\Temp\ViralGenomeProbes_withOverlap.fasta')    

        %Andi's adapter sequences
        adapterA_5='ATACGCCATGCTGGGTCTCC';
        adapterA_3='CGTACTTGGGAGTCGGCCAT';
        adapterB_5='AGGCCCTGGCTGCTGATATG';
        adapterB_3='GACCTTTTGGGACAGCGGTG';

        [ res ] = makeViralProbesAlternatingAdapters( sequenceFileName, outFileName, ...
            adapterA_5, adapterA_3, adapterB_5, adapterB_3, overlap, fragLength, rc);
        
        
        This will make 100 bp tiling fragments along the genomes in the
        C:\Temp\ViralGenome.fasta file that start every 50 bases.
        
        NEXT STEP 1
        This does not remove duplicates.  From here, to remove duplicates
        you have to run Tim Fennel's RemoveSimilarSequences.jar
        
        IF you want to use the offset you need to run this and then add the
        primers later.  It considers the primers part of the probe so it
        doesn't work for that yet, so you have to modify the code or send
        it blank adapters and then put the adapters on later.  You can
        probably do this with a simple unix script.
        
        example:
        java -jar ~tfennell/bin/RemoveSimilarSequences.jar OFFSET_ALLOWED=20 MISMATCHES_ALLOWED=4 I= C:\Temp\ViralGenomeProbes_withOverlap.fasta O=Probes_unique4.fasta
        
        NEXT STEP 2
        Then Tarjei doesn't want the header sequences.  He just wants the
        probes so you can remove the > lines with a grep
        e.g. grep -v '>' Probes_unique4.fasta > justSeqs.txt

        
        
        
      %}


    
    display('reading')
    outFileName
    
    
    [ seqNames, seqs]= ReadFasta(sequenceFileName);
    fid=fopen( outFileName, 'w');
    
    LSeqs=length(seqs)
    
    fragLength

    for i=1:length(seqs)

        thisSeq=seqs{i};
        fragCtr=1;    
        
      

        for ctr=1:fragLength:length(thisSeq)-fragLength    
           
            %The A adapter
            thisFrag=thisSeq(ctr:ctr+fragLength-1);

            %Replace the degenerate bases with N
            thisFrag=regexprep(thisFrag, '[YRWSMK]', 'N');

            [matchstart,~,~,~,~,~,~]=regexp(thisFrag, 'N{2,}')     ;
            
            
            if isempty(matchstart)        
                probeSeq=strcat(adapterA_5, thisFrag, adapterA_3);
                fprintf(fid,'%s%s%s%s\n', '>', seqNames{i}, '_Forward_', num2str(fragCtr));  
                fprintf(fid,'%s\n', probeSeq);  
                if rc==1
                  fprintf(fid,'%s%s%s%s\n', '>', seqNames{i}, '_RC_', num2str(fragCtr));  
                  fprintf(fid,'%s\n',seqrcomplement(probeSeq));
                end
            end

            fragCtr=fragCtr+1;

            %The B adapter
            if length(thisSeq)>=ctr+overlap+fragLength-1

                thisFrag=thisSeq(ctr+overlap:ctr+overlap+fragLength-1);

                %Replace the degenerate bases with N
                thisFrag=regexprep(thisFrag, '[YRWSMK]', 'N');
                [matchstart,~,~,~,~,~,~]=regexp(thisFrag, 'N{1,}');
                
                if isempty(matchstart) 
                    probeSeq=strcat(adapterB_5, thisFrag, adapterB_3);                
                    fprintf(fid,'%s%s%s%s\n', '>', seqNames{i}, '_Forward_', num2str(fragCtr));  
                    fprintf(fid,'%s\n', probeSeq);  
                    
                    if rc==1
                      fprintf(fid,'%s%s%s%s\n', '>', seqNames{i}, '_RC_', num2str(fragCtr));  
                      fprintf(fid,'%s\n',seqrcomplement(probeSeq)); 
                    end
                end

                 fragCtr=fragCtr+1;
            end

        end

    end


    %Add in reads adjacent to strings of Ns

    ctr=1;
    for i=1:length(seqs)
        s=seqs{i};
        s=regexprep(s, '[YRWSMK]', 'N');
        [matchstart,matchend,~,~,~,~,splitstring]=regexp(s, 'N{2,}');
        for j=1:length(matchstart)
           %make a left probe
           if matchstart(j)>fragLength           
               possibleProbe=s(matchstart(j)-fragLength:matchstart(j)-1);
               if size(regexp(possibleProbe, 'N{2,}'))==0           
                   ctr=ctr+1;
                   x=strcat( adapterA_3, possibleProbe, adapterA_5);
                   probe{1}=x;               
                   x=strSplit(seqNames{i}, ' ');
                   y=strcat(x(1), '_bait_Left_of_Ns_at_', num2str(matchstart(j)));
                   probeName{1}=y{1};           
                   WriteFasta( outFileName,probeName, probe, 1,2000)       

                   ctr=ctr+1;
                    
                   if(rc==1)
                       x= seqrcomplement(strcat( adapterA_5, possibleProbe, adapterA_5));
                       probe{1}=x;               
                       x=strSplit(seqNames{i}, ' ');
                       y=strcat(x(1), '_bait_Left_of_Ns_at_', num2str(matchstart(j)), '_RC');
                       probeName{1}=y{1};           
                       WriteFasta( outFileName,probeName, probe, 1,2000)
                   end

               end
           end

           if length(s)-matchend(j)>fragLength   
               display('running right')
               possibleProbe=s(matchend(j)+1:matchend(j)+fragLength);
               if size(regexp(possibleProbe, 'N{2,}'))==0     
                   display(' printing right')
                   ctr=ctr+1;               
                   x=strcat( adapterA_5, possibleProbe, adapterA_5);
                   probe{1}=x;       
                   x=strSplit(seqNames{i}, ' ');
                   y=strcat(x(1), '_bait_Right_of_Ns_at_', num2str(matchstart(j)));
                   probeName{1}=y{1};   
                   WriteFasta( outFileName,probeName, probe, 1,2000)
                    
                   if rc==1
                       x= seqrcomplement(strcat( adapterA_5, possibleProbe, adapterA_3));
                       probe{1}=x;       
                       x=strSplit(seqNames{i}, ' ');
                       y=strcat(x(1), '_bait_Right_of_Ns_at_', num2str(matchstart(j)), '_RC');
                       probeName{1}=y{1};   
                       WriteFasta( outFileName,probeName, probe, 1,2000)
                   end

               end
           end

        end

    end

    res=0;

end



