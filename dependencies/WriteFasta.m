function [ seqNames, seqs] = WriteFasta( outFileName, seqNames, seqs, append, lineLength)

    
    %Write all of the sequences and seqNames that are there
    %with 80 chars per line
    if append==1
        fid=fopen( outFileName, 'a');
    else
        fid=fopen( outFileName, 'w');
    end 
    
    for i=1:length(seqNames)
      
       fprintf(fid,'>%s\n',seqNames{i});
       
       seqStr=seqs{i};
       while isempty(seqStr)==0
           L=length(seqStr);
           cut=min([L lineLength]);
           try
            fprintf(fid,'%s\n', seqStr(1:cut));    
           catch
               display('Cannot print seq')
           end
           seqStr=seqStr(cut+1:L);
       end
    end
 
    fclose(fid);        

end

