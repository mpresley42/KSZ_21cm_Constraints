! Takes in a binary data cube from 21cmFAST and outputs a single slice
! in ASCII format for visualization
! Written by Adrian Liu, 14th July 2010.
      program extractSpatialSliceAv
      implicit none
      integer boxlength
      character*60 boxfname,avTbFname

      open(2,file='args.dat',status='old')
      read(2,*,end=777,err=777) boxfname,boxlength,avTbFname
      close(2)

      print *, boxfname,boxlength,avTbFname
      call go(boxfname,boxlength,avTbFname)
      
      stop
 777  call usage
      end
      
      subroutine usage
      print *, 'Averages over a single precision binary cube'
      print *, 'USAGE: call extractSpatialAv.x <box fname> <box length>'
      print *, '<spatial average fname>'
      print *, 'EXAMPLE: call extractSpatialAv.x box.dat 300 avTb.dat'
      return
      end 
      
      subroutine go(boxfname,boxlength,avTbFname)
      implicit none
      integer i,j,k,boxlength
      real*4 field(boxlength,boxlength,boxlength)
      real*8 spatialAv
      character*60 boxfname,avTbFname

      print *, boxfname,boxlength,avTbFname
      
! Initialize the field array
      do i=1,boxlength
         do j=1,boxlength
            do k=1,boxlength
               field(i,j,k)=0.0
            enddo
         enddo
      enddo

! Read in the data cube
      open(72,file=boxfname,form='unformatted',access='direct',recl=4)
      do k=1,boxlength
         do j=1,boxlength
            do i=1,boxlength
               read(72,rec=(k-1)*boxlength*boxlength+(j-1)*boxlength+i) field(i,j,k)
            enddo
         enddo
      enddo
      close(72)

! Compute the spatial average for the chosen slice
      spatialAv=0.0
      do i=1,boxlength
         do j=1,boxlength
            do k=1,boxlength
               spatialAv=spatialAv+field(i,j,k)
            enddo
         enddo
      enddo
      spatialAv=spatialAv/(boxlength*boxlength*boxlength)
      print *, spatialAv
      print *, field(45,45,45)
      print *, avTbFname
      open(28,file=avTbFname)
      write(28,*) spatialAv
      close(28)

      return
      end
