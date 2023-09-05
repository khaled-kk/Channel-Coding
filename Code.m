clear all
close all
clc
obj=VideoReader('highway.avi'); 
a=read(obj);

frames=get(obj,'NumberOfFrames');
for i=1:frames             %30 frames
    I(i).cdata=a(:,:,:,i); %Red, green, blue for each frame
                           %cdata >>Allocating the positions of the colours
end                        %I>> Image, (.cdata: data of the image)
               

s=size(I(1).cdata); 
mov(1:frames) =struct('cdata', zeros(s(1),s(2), 3, 'uint8'),'colormap', []); 
%generating a new video with the same size of the original one


for i=1:frames 
    
 FinalArray = [];     %concatinated decoded sequence
   
 R=I(i).cdata(:,:,1); 
 G=I(i).cdata(:,:,2);
 B=I(i).cdata(:,:,3); %extract the colour of each colour in order
 

 Rdouble = double(R);
 Gdouble = double(G);
 Bdouble = double(B); % 1 row >> 144 x 176
 
 Rbin = de2bi(Rdouble);
 Gbin = de2bi(Gdouble);
 Bbin = de2bi(Bdouble); %converting from Decimal to Binary, 144 x 176 rows and 8 coloumns
 
 Rbin1 = reshape(Rbin, 1 ,202752); 
 Gbin1 = reshape(Gbin, 1 ,202752);
 Bbin1 = reshape(Bbin, 1 ,202752);
 
 WholeFrame= [Rbin1 Gbin1 Bbin1]; %putting the colours altogether in 1 row

 
 trellis = poly2trellis(7,[171 133]);
 %ConstraintLength ,CodeGenerator
 %Convert convolutional code polynomials to trellis description
 
    Puncturing_rule1 = [1 1 1 0 1 0 1 0 0 1 1 0 1 0 1 0]; %8/9
    Puncturing_rule2 = [1 1 1 0 1 0 1 0 1 1 1 0 1 0 1 0]; %4/5
    Puncturing_rule3 = [1 1 1 0 1 1 1 0 1 1 1 0 1 1 1 0]; %2/3
    Puncturing_rule4 = [1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 0]; %4/7
    Puncturing_rule5 = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]; %1/2

    for k=1:1024:length(WholeFrame) %each time increments by 1024(size of the packet) to loop on the packets
        
          %decoded sequence are the whole frame bits
          
        Encoded = convenc(WholeFrame(k:k+1023), trellis , Puncturing_rule1);                                                                     
                                                                             
        Errored = bsc(Encoded,0.01 ); %calculates the errored from the encoded by probabilty of error 0.01
          
        DecodedSequence=vitdec(Errored,trellis, 35 ,'trunc','hard',Puncturing_rule1); %35>> depth and +ve#
                                               
        if(DecodedSequence==WholeFrame(k:k+1023)) %checking if the decoded is equal to the original to make sure the puncturing used is valid
            
             FinalArray = cat(2,FinalArray,DecodedSequence);  %concatinated decoded sequence
            
        continue  %goes directly to the next iteration and ends the if we're in >> +1024
        
        end    
   
        Encoded  = convenc(WholeFrame(k:k+1023), trellis , Puncturing_rule2); %4/5
          
        Errored  = bsc(Encoded,0.01 );
          
        DecodedSequence=vitdec(Errored,trellis, 35 ,'trunc','hard',Puncturing_rule2);
        
          
        if(DecodedSequence==WholeFrame(k:k+1023))
            
            FinalArray = cat(2,FinalArray,DecodedSequence);
            
        continue 
        
        end    
           
        
        Encoded  = convenc(WholeFrame(k:k+1023), trellis , Puncturing_rule3); 
          
        Errored  = bsc(Encoded,0.01 );
          
        DecodedSequence=vitdec(Errored,trellis, 35 ,'trunc','hard',Puncturing_rule3);
        
        
        
        if(DecodedSequence==WholeFrame(k:k+1023))
            
         FinalArray = cat(2,FinalArray,DecodedSequence);  
            
        continue   
        
        end    
            
        
        Encoded  = convenc(WholeFrame(k:k+1023), trellis , Puncturing_rule4); %4/7
          
        Errored  = bsc(Encoded,0.01 );
          
        DecodedSequence=vitdec(Errored,trellis, 35 ,'trunc','hard',Puncturing_rule4);
        
        if(DecodedSequence==WholeFrame(k:k+1023))
            
            FinalArray = cat(2,FinalArray,DecodedSequence);
            
        continue   
        end    
        
        Encoded  = convenc(WholeFrame(k:k+1023), trellis , Puncturing_rule5);
          
        Errored  = bsc(Encoded,0.01 );
          
        DecodedSequence=vitdec(Errored,trellis, 35 ,'trunc','hard',Puncturing_rule5); 
        
        FinalArray = cat(2,FinalArray,DecodedSequence);
        
        continue
 
    end
    
        RedOne = FinalArray(1:202752); 
        GreenOne = FinalArray(202753:405504);
        BlueOne = FinalArray (405505:608256); %each range resembles the numbrt of coloumns per colour
        
        RedTwo = reshape(RedOne,25344,8);     %25344 = 144 x 176
        GreenTwo = reshape(GreenOne,25344,8);
        BlueTwo = reshape (BlueOne,25344,8); %middle=no.of rows, last=no.of coloumns
        
        RedDouble= bi2de(RedTwo);
        GreenDouble= bi2de(GreenTwo);
        BlueDouble= bi2de(BlueTwo);
        
        Reduint8 = uint8(RedDouble);
        Greenuint8 = uint8(GreenDouble);
        Blueuint8 = uint8(BlueDouble); %change to unsigned 8 bit integar
        
        RedThree = reshape (Reduint8, 144,176);
        GreenThree = reshape (Greenuint8, 144,176);
        BlueThree = reshape (Blueuint8, 144,176); 
        
        Rbin2 = RedThree;
        Gbin2 = GreenThree;
        Bbin2 = BlueThree; 
        
mov(1,i).cdata(:,:,1) = Rbin2; %constructing a frame specifying the position of the colour in the cdata format
mov(1,i).cdata(:,:,2) = Gbin2;
mov(1,i).cdata(:,:,3) = Bbin2;

end 

obj_new = VideoWriter('NewVideoFile.avi', 'motion JPEG AVI');
open(obj_new)
writeVideo(obj_new,mov)
close(obj_new)
implay('NewVideoFile.avi')  









