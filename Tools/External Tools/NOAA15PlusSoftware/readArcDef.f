      program readArcDef(arcfile,outfile)
      
      ! For reading an archiveDefaultsSS file to check date.
      ! The program takes two arguments, specifying the input
      ! filename and the output filename (including paths.)  
      ! The output file is written in the current working directory 
      ! unless a path is included as part of the filename.
      
      integer*2 j
      integer*4 rec,eof,hr,min,sec
      character*80 arcfile,outfile
      
      ! Common block /rec/ variables:
      integer*4 cSum,ihd(6,4)
      integer*2 cSumFlag,major,status(10),qual(16),minor(16),mdf(40,16)
      real*4 analog(17),head(27,4),ssLoc(2,16),mep0(9,16),mep90(9,16),
     +    mepOmni(4,16),ted0(8,16),ted30(8,16),ted0s(8,4),ted30s(8,4),
     +    tedback(2,4),tedfx(7,16)
      
      common /rec/cSum,ihd,cSumFlag,major,status,qual,minor,
     +    mdf,analog,head,ssLoc,mep0,mep90,mepOmni,ted0,ted30,
     +    ted0s,ted30s,tedback,tedfx
      
      eof = 0
      open(unit=7,access='direct',recl=2544,file = arcfile,status='old')
      open(unit=11,access='sequential',file = outfile,
     +     status='unknown',form='formatted')
      
      do rec = 1,1  
      call readArcSub(7,rec,eof)
      if(eof.eq.1) goto 100
      do j=1,4 
          hr = ihd(4,j)/3600000
          min = (ihd(4,j) - hr*3600000)/60000
          sec = (ihd(4,j) - hr*3600000 - min*60000)/1000              
          write(11,200),ihd(2,j),ihd(3,j),hr,min,sec,
     +                  head(2,j),head(3,j),ihd(6,j)              
      enddo
      enddo
  100 close(7)
      close(11)
  200 format(i4,1x,i3,3(1x,i2),1x,f8.3,1x,f8.3,1x,i7)      
      end          
