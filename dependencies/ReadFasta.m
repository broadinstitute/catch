function [ seqNames, seqs] = ReadFasta( inFileName)

    %Not mem efficient, see if it works
    fid=fopen( inFileName, 'r');
      
    seqCtr=0;

    %Read file line by line and get sequenceNames and sequences
    oldSeq='';
    tline = fgetl(fid);    
    while ischar(tline)            
        if isempty(tline)==0
            if strcmp(tline(1),'>')==1
                if seqCtr>0
                    seqs{seqCtr}=oldSeq;
                end
                seqCtr=seqCtr+1;
                seqNames{seqCtr}=tline(2:length(tline));                
                oldSeq='';
            else
                oldSeq=strcat(oldSeq, tline);
            end
        end 
        tline = fgetl(fid);    
    end

    seqs{seqCtr}=oldSeq;
    
    fclose(fid);        

end

