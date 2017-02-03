! Takes in two binary data cubes from 21cmFAST:
! i) Density cubes (which are no longer zero-padded in latest 21cmFAST versions)
! ii) Neutral fraction cubes (which are not)
! Pick a slice, and perform a spatial average of the
! quantity x_e (1+delta), where x_e = (1-x_H) is the
! ionized fraction.  Also outputs the slice
! Written by Adrian Liu, 29th May 2014.
      program compute_weightedIonizedFrac_fullBox
      implicit none
      integer boxlength
      character*60 densityFname,neutralFname,weightedxeFname

      open(2,file='args.dat',status='old')
      read(2,*,end=777,err=777) densityFname,neutralFname,boxlength,weightedxeFname
      close(2)

!     print *, boxlength
      call go(densityFname,neutralFname,boxlength,weightedxeFname)
      
      stop
 777  call usage
      end
      
      subroutine usage
      print *, 'Outputs x_e (1+delta), spatially averaged over the whole box'
      print *, 'USAGE: call compute_weightedIonizedFrac.x <density box fname> <neutral frac fname>'
      print *, '<box length> <spatial average fname>'
      print *, 'EXAMPLE: call compute_weightedIonizedFrac.x density_box.dat ion_box.dat'
      print *, '200 av_weightedxe.dat'
      return
      end 
      
      subroutine go(densityFname,neutralFname,boxlength,weightedxeFname)
      implicit none
      integer i,j,k,boxlength,adjLen
      real*4 densityField(boxlength,boxlength,boxlength)
      real*4 ionizedFracField(boxlength,boxlength,boxlength)
      real*8 spatialAv
      character*60 densityFname,weightedxeFname,neutralFname
      
! Initialize the field arrays
      do i=1,boxlength
         do j=1,boxlength
            do k=1,boxlength
               densityField(i,j,k)=0.0
               ionizedFracField(i,j,k)=0.0
            enddo
         enddo
      enddo

! Read in the data cubes
      open(72,file=densityFname,form='unformatted',access='direct',recl=4)
      open(73,file=neutralFname,form='unformatted',access='direct',recl=4)
      do i=1,boxlength
         do j=1,boxlength
            do k=1,boxlength
               adjLen=boxlength !+2
               read(72,rec=(i-1)*adjLen*boxlength+(j-1)*adjLen+k) densityField(i,j,k)
               read(73,rec=(i-1)*boxlength*boxlength+(j-1)*boxlength+k) ionizedFracField(i,j,k)
               ionizedFracField(i,j,k)=1.0-ionizedFracField(i,j,k) ! The input is neutral frac
            enddo
         enddo
      enddo
      close(72)
      close(73)

! Compute the spatial average for the chosen slice
      spatialAv=0.0
      do i=1,boxlength
         do j=1,boxlength
            do k=1,boxlength
               if (densityField(i,j,k) .le. -1.0) then
                  print *, "Warning: negative density detected"
               endif
               spatialAv=spatialAv+(1.0+densityField(i,j,k))*ionizedFracField(i,j,k)
            enddo
         enddo
      enddo
      spatialAv=spatialAv/(boxlength*boxlength*boxlength)
      open(28,file=weightedxeFname)
      write(28,*) spatialAv
      close(28)

      return
      end
