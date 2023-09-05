clear all
close all
clc
obj=VideoReader('foreman.avi'); 
a=read(obj);

frames=get(obj,'NumberOfFrames');
for i=1:frames
    I(i).cdata=a(:,:,:,i);
end

s=size(I(1).cdata);
mov(1:frames) =struct('cdata', zeros(s(1),s(2), 3, 'uint8'),'colormap', []);

error=0.0001:0.01:0.2;
TotalNeededData=zeros(1,max(size(error)));
TotalSentData=zeros(1,max(size(error)));
BitsInError=zeros(1,max(size(error)));
ProbabiltyOfError=zeros(1,max(size(error)));
z=1;
counter=1;


for i=1:max(size(error))
 
    x=length(i);
    TotalErroredBits=0;
    FinalArray = [];
 R=I(1).cdata(:,:,1);
 G=I(1).cdata(:,:,2);
 B=I(1).cdata(:,:,3);
%  R1=reshape(R, 1, []);
%  G1=reshape(G, 1, []);
%  B1=reshape(B, 1, []);
 Rdouble = double(R);
 Gdouble = double(G);
 Bdouble = double(B);
 Rbin = de2bi(Rdouble);
 Gbin = de2bi(Gdouble);
 Bbin = de2bi(Bdouble);
 Rbin1 = reshape(Rbin, 1 ,202752); %middle=no.of rows
 Gbin1 = reshape(Gbin, 1 ,202752);
 Bbin1 = reshape(Bbin, 1 ,202752);
 WholeFrame= [Rbin1 Gbin1 Bbin1];
 
 trellis = poly2trellis(7,[171 133]);


    Puncturing_rule1 = [1 1 1 0 1 0 1 0 0 1 1 0 1 0 1 0]; %for code rate 8/9
    Puncturing_rule2 = [1 1 1 0 1 0 1 0 1 1 1 0 1 0 1 0]; %for code rate 4/5
    Puncturing_rule3 = [1 1 1 0 1 1 1 0 1 1 1 0 1 1 1 0]; %for code rate 2/3
    Puncturing_rule4 = [1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 0]; %for code rate 4/7
    Puncturing_rule5 = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]; %for code rate 1/2
    
    
    ProbabiltyOfError(z) = i;
    z=z+1;

    for k=1:1024:length(WholeFrame)
       
        disp(k);
        Encoded = convenc(WholeFrame(k:k+1023), trellis , Puncturing_rule1);
          
        Errored = bsc(Encoded,error(i) );
          
        DecodedSequence=vitdec(Errored,trellis, 35 ,'trunc','hard',Puncturing_rule1);
        
        if(DecodedSequence==WholeFrame(k:k+1023))
             
            
             TotalSentData(i) = TotalSentData(i)+1152;
             TotalNeededData(i)=TotalNeededData(i)+1024;
             FinalArray = cat(2,FinalArray,DecodedSequence);
             counter = counter+1;
            
        continue  
        
        end    
   
        Encoded  = convenc(WholeFrame(k:k+1023), trellis , Puncturing_rule2);
          
        Errored  = bsc(Encoded,error(i) );
          
        DecodedSequence=vitdec(Errored,trellis, 35 ,'trunc','hard',Puncturing_rule2);
        
          
        if(DecodedSequence==WholeFrame(k:k+1023))
            
              TotalSentData(i) = TotalSentData(i)+1280;
              TotalNeededData(i)=TotalNeededData(i)+1024;
              FinalArray = cat(2,FinalArray,DecodedSequence);
              counter = counter+1;
            
        continue 
        
        end    
           
        
        Encoded  = convenc(WholeFrame(k:k+1023), trellis , Puncturing_rule3);
          
        Errored  = bsc(Encoded,error(i) );
          
        DecodedSequence=vitdec(Errored,trellis, 35 ,'trunc','hard',Puncturing_rule3);
        
        
        if(DecodedSequence==WholeFrame(k:k+1023))
            
             TotalSentData(i) = TotalSentData(i)+1536;
             TotalNeededData(i)=TotalNeededData(i)+1024;
             FinalArray = cat(2,FinalArray,DecodedSequence);
             counter = counter+1;
       
        continue   
        
        end    
            
        
        Encoded  = convenc(WholeFrame(k:k+1023), trellis , Puncturing_rule4);
          
        Errored  = bsc(Encoded,error(i) );
          
        DecodedSequence=vitdec(Errored,trellis, 35 ,'trunc','hard',Puncturing_rule4);
        
        if(DecodedSequence==WholeFrame(k:k+1023))
            
               TotalSentData(i) = TotalSentData(i)+1792;
               TotalNeededData(i)=TotalNeededData(i)+1024;
               FinalArray = cat(2,FinalArray,DecodedSequence);
               counter = counter+1;
            
           
            
        continue   
        end    
        
         Encoded  = convenc(WholeFrame(k:k+1023), trellis , Puncturing_rule5);
          
         Errored  = bsc(Encoded,error(i) );
          
         DecodedSequence=vitdec(Errored,trellis, 35 ,'trunc','hard',Puncturing_rule5);
        
         FinalArray = cat(2,FinalArray,DecodedSequence);
        
          TotalSentData(i) = TotalSentData(i)+2048;
          TotalNeededData(i)=TotalNeededData(i)+ 1024;
          counter = counter+1;
       
        
        continue
 
    end
         TotalErroredBits= TotalErroredBits + sum(xor(WholeFrame,FinalArray));
         BitsInError(i) =  TotalErroredBits;
         TotalErroredBits =0;
         
       
end 
 
 BER= BitsInError./608256;
 plot(error, BER)
 title('BER vs p error values for incremental redundancy');
 xlabel("Error");
 ylabel("BER");
 