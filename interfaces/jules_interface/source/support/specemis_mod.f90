!----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!----------------------------------------------------------------------------
!> @brief Calculates a spectrally varying emissivity for a given set of bands.
!>        A Planck weighting is applied. Excluded bands are taken care of.
!>        In addition, a grey emissivity across the longwave is calculated.
!>        Data for the spectrally varying emissivities are available
!>        currently for sea and desert tiles from the Feldman paper:
!>        Daniel R. Feldman, William D. Collins, Robert Pincus, Xianglei Huang,
!>        and Xiuhong Chen, Far-infrared surface emissivity and climate,
!>        PNAS November 18, 2014 111 (46) 16297-16302;
!>        https://doi.org/10.1073/pnas.1413640111.
!>        and alternatively for sea using the IREMIS model
!>        Roger Saunders, James Hocking, Emma Turner, Peter Rayer, David Rundle,
!>        Pascal Brunel, Jerome Vidot, Pascale Roquet, Marco Matricardi,
!>        Alan Geer, Niels Bormann, and Cristina Lupu,
!>        An update on the RTTOV fast radiative transfer model
!>        (currently at version 12), Geosci. Model Dev., 11, 2717-2737, 2018;
!>        https://doi.org/10.5194/gmd-11-2717-2018.
!>        The emissivity data from the paper and the model are contained in
!>        lookup tables in this routine. The data from the paper were
!>        obtained by scanning in the relevant figures in the paper.

module specemis_mod

  implicit none

  contains

! @param[in]  surf_type          String to identify the surface emissivity spectrum (lookup table)
! @param[in[  n_bands            Number of spectral bands
! @param[in]  wavelen_low        Low wavelength limits of the spectral bands (units: metre)
! @param[in]  wavelen_high       High wavelength limits of the spectral bands (units: metre)
! @param[in[  tile_temp          Tile temperature (required in Planck weighting, units: Kelvin)
! @param[in]  n_band_exclude     Number of excluded bands for the spectral bands
! @param[in]  index_band_exclude Indices of the excluded bands for the spectral bands
! @param[out] emis               Emissivities for each of the spectral bands (dimensionless)
! @param[out] grey_emis          Grey emissivity for the longwave region (dimensionless)
  subroutine specemis(surf_type, n_bands, wavelen_low, wavelen_high, &
             tile_temp, n_band_exclude, index_band_exclude, emis,grey_emis)

  use constants_mod, only : r_def, i_def

  use planck_mod, only : planck

  implicit none

  ! String to identify the surface emissivity spectrum (lookup table)
  character(len=*), intent(in) :: surf_type
  ! Number of spectral bands
  integer(i_def), intent(in) :: n_bands
  ! Low wavelength limits of the spectral bands (units: metre)
  real(r_def), intent(in) :: wavelen_low(n_bands)
  ! High wavelength limits of the spectral bands (units: metre)
  real(r_def), intent(in) :: wavelen_high(n_bands)
  ! Tile temperature (required in Planck weighting, units: Kelvin)
  real(r_def), intent(in) :: tile_temp
  ! Number of excluded bands for the spectral bands
  integer(i_def), intent(in) :: n_band_exclude(n_bands)
  ! Indices of the excluded bands for the spectral bands
  integer(i_def), intent(in) :: index_band_exclude(:,:)
  ! Emissivities for each of the spectral bands (dimensionless)
  real(r_def), intent(out) :: emis(n_bands)
  ! Grey emissivity for the longwave region (dimensionless)
  real(r_def), intent(out) :: grey_emis

  ! Size of lookup table (0-3000 wavenumbers at 10 wavenumber spacing)
  ! Units of wavenumbers are inverse cm
  integer(i_def), parameter :: n_val=301
  ! Planck function at lookup table wavelenghts
  real(r_def) :: planck_wavelentbl(n_val)

  ! Planck weighted contribution of each lookup table position to each band
  real(r_def) :: weight(n_val,n_bands)
  ! Total weight for each band
  real(r_def) :: total_weight(n_bands)
  ! Planck weighted contribution to the grey emissivity from each lookup table position
  real(r_def) :: grey_weight(n_val)
  ! Total weight for grey emissivity
  real(r_def) :: grey_total_weight

  ! Local variables for indexing
  integer(i_def) :: i
  integer(i_def) :: j
  integer(i_def) :: k


  ! Emissivity lookup tables to hold the available literature and model data
  ! The literature data were obtained by scanning the relevant figures
  ! 0-3000 wavenumbers at 10 wavenumber spacing
  ! Units of wavenumbers are inverse cm
  ! Units of wavelengths are metre
  ! Emissivity is dimensionless

  ! Representative (centre) wavenumbers for lookup table
  real(r_def), parameter :: wavenumtbl(1:n_val)      = &
      [ 3.0_r_def, ( 10.0_r_def * ( i - 1 ), i = 2, n_val ) ]
  ! Lower wavenumber limits for lookup table
  real(r_def), parameter :: wavenumtbl_low(1:n_val)  = &
      [ 0.99_r_def, ( wavenumtbl(i) - 5.0_r_def, i = 2, n_val ) ]
  ! High wavenumber limits for lookup table
  real(r_def), parameter :: wavenumtbl_high(1:n_val) = &
      [ 5.0_r_def, ( wavenumtbl(i) + 5.0_r_def, i = 2, n_val ) ]

  ! Lower wavelengths limits for lookup table (units: metre)
  real(r_def), parameter :: wavelentbl_low(1:n_val)  = &
      [ ( 1.0e-2_r_def / wavenumtbl_high(i), i = 1, n_val ) ]
  ! Higher wavelengths limts for lookup table (units: metre)
  real(r_def), parameter :: wavelentbl_high(1:n_val) = &
      [ ( 1.0e-2_r_def / wavenumtbl_low(i), i = 1, n_val ) ]
  ! Representative wavelengths for lookup table (units: metre)
  real(r_def), parameter :: wavelentbl(1:n_val)      = &
      [ (0.5_r_def * ( wavelentbl_low(i) + wavelentbl_high(i) ), i = 1, n_val ) ]

  ! Emissivity spectrum for desert from Feldman
  real(r_def), parameter :: emis_desert_feldman(1:n_val) = &
    (/ 0.911264_r_def, 0.911264_r_def, 0.911264_r_def, 0.911264_r_def, 0.911104_r_def, &
       0.910943_r_def, 0.910943_r_def, 0.910943_r_def, 0.910621_r_def, 0.910621_r_def, &
       0.910621_r_def, 0.911586_r_def, 0.913515_r_def, 0.914801_r_def, 0.913997_r_def, &
       0.913676_r_def, 0.912872_r_def, 0.913354_r_def, 0.914158_r_def, 0.915444_r_def, &
       0.917212_r_def, 0.915605_r_def, 0.913837_r_def, 0.914480_r_def, 0.914640_r_def, &
       0.913033_r_def, 0.911425_r_def, 0.909175_r_def, 0.910943_r_def, 0.912711_r_def, &
       0.913354_r_def, 0.911747_r_def, 0.908049_r_def, 0.906603_r_def, 0.906120_r_def, &
       0.909335_r_def, 0.910621_r_def, 0.908692_r_def, 0.906603_r_def, 0.904191_r_def, &
       0.897440_r_def, 0.891331_r_def, 0.884901_r_def, 0.880721_r_def, 0.881043_r_def, &
       0.872845_r_def, 0.861431_r_def, 0.851947_r_def, 0.869790_r_def, 0.886026_r_def, &
       0.893421_r_def, 0.890206_r_def, 0.884740_r_def, 0.877346_r_def, 0.877346_r_def, &
       0.885062_r_def, 0.892617_r_def, 0.904674_r_def, 0.914801_r_def, 0.924125_r_def, &
       0.930716_r_def, 0.937949_r_def, 0.943576_r_def, 0.950810_r_def, 0.954828_r_def, &
       0.956973_r_def, 0.959330_r_def, 0.961580_r_def, 0.962545_r_def, 0.964795_r_def, &
       0.966563_r_def, 0.967000_r_def, 0.967367_r_def, 0.964474_r_def, 0.963348_r_def, &
       0.962062_r_def, 0.961580_r_def, 0.960294_r_def, 0.958847_r_def, 0.957240_r_def, &
       0.956597_r_def, 0.955632_r_def, 0.954507_r_def, 0.953542_r_def, 0.951613_r_def, &
       0.950167_r_def, 0.947755_r_def, 0.945826_r_def, 0.943254_r_def, 0.940682_r_def, &
       0.936181_r_def, 0.932323_r_def, 0.932805_r_def, 0.934252_r_def, 0.932645_r_def, &
       0.931198_r_def, 0.927500_r_def, 0.920588_r_def, 0.914480_r_def, 0.905156_r_def, &
       0.896797_r_def, 0.885383_r_def, 0.880078_r_def, 0.871237_r_def, 0.860627_r_def, &
       0.858377_r_def, 0.862556_r_def, 0.872041_r_def, 0.886509_r_def, 0.896797_r_def, &
       0.909175_r_def, 0.917212_r_def, 0.919945_r_def, 0.918498_r_def, 0.920588_r_def, &
       0.923803_r_def, 0.929430_r_def, 0.935699_r_def, 0.946469_r_def, 0.954185_r_def, &
       0.964956_r_def, 0.972190_r_def, 0.979745_r_def, 0.984889_r_def, 0.986818_r_def, &
       0.987461_r_def, 0.986979_r_def, 0.986175_r_def, 0.985050_r_def, 0.984085_r_def, &
       0.982799_r_def, 0.982317_r_def, 0.981500_r_def, 0.980500_r_def, 0.979263_r_def, &
       0.979102_r_def, 0.978459_r_def, 0.977977_r_def, 0.977334_r_def, 0.977012_r_def, &
       0.976369_r_def, 0.975887_r_def, 0.975566_r_def, 0.975244_r_def, 0.974923_r_def, &
       0.974601_r_def, 0.974440_r_def, 0.974280_r_def, 0.974119_r_def, 0.973958_r_def, &
       0.973797_r_def, 0.973476_r_def, 0.973315_r_def, 0.972994_r_def, 0.972672_r_def, &
       0.972351_r_def, 0.972351_r_def, 0.972029_r_def, 0.971868_r_def, 0.971708_r_def, &
       0.971547_r_def, 0.971547_r_def, 0.971386_r_def, 0.971386_r_def, 0.971386_r_def, &
       0.971547_r_def, 0.971547_r_def, 0.971386_r_def, 0.971386_r_def, 0.971386_r_def, &
       0.971064_r_def, 0.971064_r_def, 0.970904_r_def, 0.970743_r_def, 0.970743_r_def, &
       0.970421_r_def, 0.970421_r_def, 0.970100_r_def, 0.970100_r_def, 0.969939_r_def, &
       0.969778_r_def, 0.969778_r_def, 0.969778_r_def, 0.969457_r_def, 0.969457_r_def, &
       0.969457_r_def, 0.969296_r_def, 0.969135_r_def, 0.969135_r_def, 0.968975_r_def, &
       0.968975_r_def, 0.968814_r_def, 0.968814_r_def, 0.968814_r_def, 0.968492_r_def, &
       0.968492_r_def, 0.968492_r_def, 0.968332_r_def, 0.968492_r_def, 0.968171_r_def, &
       0.968171_r_def, 0.968171_r_def, 0.968010_r_def, 0.967849_r_def, 0.967849_r_def, &
       0.968171_r_def, 0.967849_r_def, 0.967689_r_def, 0.967689_r_def, 0.967689_r_def, &
       0.967689_r_def, 0.967689_r_def, 0.967689_r_def, 0.967689_r_def, 0.967689_r_def, &
       0.967528_r_def, 0.967528_r_def, 0.967206_r_def, 0.967206_r_def, 0.967206_r_def, &
       0.967046_r_def, 0.966885_r_def, 0.966885_r_def, 0.967046_r_def, 0.966885_r_def, &
       0.967046_r_def, 0.966724_r_def, 0.966724_r_def, 0.966724_r_def, 0.966724_r_def, &
       0.966724_r_def, 0.966403_r_def, 0.966242_r_def, 0.966242_r_def, 0.966081_r_def, &
       0.965920_r_def, 0.965920_r_def, 0.965920_r_def, 0.965760_r_def, 0.965599_r_def, &
       0.965599_r_def, 0.965599_r_def, 0.965438_r_def, 0.965438_r_def, 0.965277_r_def, &
       0.965277_r_def, 0.965277_r_def, 0.965117_r_def, 0.965117_r_def, 0.964956_r_def, &
       0.964956_r_def, 0.964795_r_def, 0.964795_r_def, 0.964634_r_def, 0.964634_r_def, &
       0.964634_r_def, 0.964474_r_def, 0.964474_r_def, 0.964313_r_def, 0.964313_r_def, &
       0.964313_r_def, 0.964152_r_def, 0.963991_r_def, 0.963991_r_def, 0.963991_r_def, &
       0.963831_r_def, 0.963831_r_def, 0.963670_r_def, 0.963670_r_def, 0.963670_r_def, &
       0.963509_r_def, 0.963348_r_def, 0.963348_r_def, 0.963348_r_def, 0.963188_r_def, &
       0.963188_r_def, 0.963348_r_def, 0.963348_r_def, 0.963348_r_def, 0.963188_r_def, &
       0.963027_r_def, 0.963027_r_def, 0.963027_r_def, 0.962866_r_def, 0.962705_r_def, &
       0.962705_r_def, 0.962705_r_def, 0.962705_r_def, 0.962545_r_def, 0.962384_r_def, &
       0.962384_r_def, 0.962384_r_def, 0.962223_r_def, 0.962223_r_def, 0.962062_r_def, &
       0.962062_r_def, 0.962062_r_def, 0.961902_r_def, 0.961741_r_def, 0.961741_r_def, &
       0.961741_r_def /)

  ! Emissivity spectrum for sea from Feldman
  real(r_def), parameter :: emis_sea_feldman(1:n_val) = &
    (/ 0.761424_r_def, 0.774032_r_def, 0.802157_r_def, 0.806197_r_def, 0.809915_r_def, &
       0.812178_r_def, 0.814118_r_def, 0.815572_r_def, 0.816219_r_def, 0.817997_r_def, &
       0.819613_r_def, 0.823331_r_def, 0.825594_r_def, 0.828665_r_def, 0.830766_r_def, &
       0.835615_r_def, 0.838848_r_def, 0.843051_r_def, 0.849354_r_def, 0.853880_r_def, &
       0.859861_r_def, 0.864387_r_def, 0.870852_r_def, 0.874408_r_def, 0.878126_r_def, &
       0.881035_r_def, 0.882975_r_def, 0.883945_r_def, 0.885076_r_def, 0.885884_r_def, &
       0.885561_r_def, 0.885399_r_def, 0.884915_r_def, 0.884915_r_def, 0.884591_r_def, &
       0.884106_r_def, 0.883783_r_def, 0.884106_r_def, 0.884106_r_def, 0.883945_r_def, &
       0.883945_r_def, 0.884268_r_def, 0.884046_r_def, 0.884137_r_def, 0.884415_r_def, &
       0.884076_r_def, 0.884238_r_def, 0.884561_r_def, 0.884399_r_def, 0.884137_r_def, &
       0.884238_r_def, 0.884415_r_def, 0.884399_r_def, 0.884313_r_def, 0.884591_r_def, &
       0.884975_r_def, 0.885884_r_def, 0.886693_r_def, 0.887662_r_def, 0.888471_r_def, &
       0.889764_r_def, 0.890572_r_def, 0.892188_r_def, 0.893320_r_def, 0.894936_r_def, &
       0.896320_r_def, 0.897846_r_def, 0.899139_r_def, 0.900917_r_def, 0.902371_r_def, &
       0.904311_r_def, 0.905604_r_def, 0.907544_r_def, 0.909807_r_def, 0.912554_r_def, &
       0.914656_r_def, 0.917404_r_def, 0.919666_r_def, 0.922738_r_def, 0.925324_r_def, &
       0.928880_r_def, 0.931628_r_def, 0.935669_r_def, 0.938578_r_def, 0.940983_r_def, &
       0.943791_r_def, 0.947084_r_def, 0.949731_r_def, 0.952877_r_def, 0.955873_r_def, &
       0.957166_r_def, 0.958136_r_def, 0.958783_r_def, 0.958621_r_def, 0.958459_r_def, &
       0.958298_r_def, 0.957166_r_def, 0.956520_r_def, 0.955388_r_def, 0.954257_r_def, &
       0.953610_r_def, 0.952479_r_def, 0.951994_r_def, 0.950701_r_def, 0.950216_r_def, &
       0.949408_r_def, 0.948761_r_def, 0.948599_r_def, 0.948115_r_def, 0.947468_r_def, &
       0.946983_r_def, 0.946821_r_def, 0.946013_r_def, 0.945690_r_def, 0.945205_r_def, &
       0.944882_r_def, 0.944397_r_def, 0.944074_r_def, 0.943104_r_def, 0.943152_r_def, &
       0.943104_r_def, 0.942781_r_def, 0.942457_r_def, 0.942296_r_def, 0.942134_r_def, &
       0.941811_r_def, 0.941487_r_def, 0.941164_r_def, 0.940841_r_def, 0.940679_r_def, &
       0.940518_r_def, 0.940194_r_def, 0.940033_r_def, 0.939709_r_def, 0.939548_r_def, &
       0.939225_r_def, 0.939225_r_def, 0.938901_r_def, 0.938578_r_def, 0.938255_r_def, &
       0.938255_r_def, 0.937931_r_def, 0.937608_r_def, 0.937285_r_def, 0.936962_r_def, &
       0.936638_r_def, 0.936315_r_def, 0.935830_r_def, 0.935669_r_def, 0.935345_r_def, &
       0.935184_r_def, 0.934860_r_def, 0.934699_r_def, 0.934375_r_def, 0.934052_r_def, &
       0.933567_r_def, 0.932759_r_def, 0.932113_r_def, 0.931304_r_def, 0.930335_r_def, &
       0.929850_r_def, 0.929041_r_def, 0.928557_r_def, 0.930981_r_def, 0.933244_r_def, &
       0.935830_r_def, 0.938093_r_def, 0.940679_r_def, 0.941487_r_def, 0.941487_r_def, &
       0.941487_r_def, 0.941487_r_def, 0.941487_r_def, 0.941326_r_def, 0.941326_r_def, &
       0.941164_r_def, 0.941164_r_def, 0.941164_r_def, 0.941164_r_def, 0.941003_r_def, &
       0.940841_r_def, 0.941003_r_def, 0.940679_r_def, 0.940841_r_def, 0.940518_r_def, &
       0.940194_r_def, 0.940033_r_def, 0.939871_r_def, 0.939548_r_def, 0.939225_r_def, &
       0.939063_r_def, 0.938901_r_def, 0.938578_r_def, 0.938255_r_def, 0.938093_r_def, &
       0.937931_r_def, 0.937608_r_def, 0.937285_r_def, 0.937123_r_def, 0.936962_r_def, &
       0.936638_r_def, 0.936315_r_def, 0.936477_r_def, 0.936315_r_def, 0.936315_r_def, &
       0.936315_r_def, 0.936315_r_def, 0.936153_r_def, 0.935992_r_def, 0.935992_r_def, &
       0.935992_r_def, 0.935992_r_def, 0.935992_r_def, 0.935830_r_def, 0.935830_r_def, &
       0.935830_r_def, 0.935669_r_def, 0.935669_r_def, 0.935669_r_def, 0.935669_r_def, &
       0.935507_r_def, 0.935345_r_def, 0.935345_r_def, 0.935345_r_def, 0.935184_r_def, &
       0.935184_r_def, 0.935022_r_def, 0.935022_r_def, 0.934699_r_def, 0.934699_r_def, &
       0.934699_r_def, 0.934537_r_def, 0.934375_r_def, 0.934375_r_def, 0.934214_r_def, &
       0.934052_r_def, 0.933891_r_def, 0.933729_r_def, 0.933729_r_def, 0.933729_r_def, &
       0.933406_r_def, 0.933406_r_def, 0.933406_r_def, 0.933244_r_def, 0.933082_r_def, &
       0.933082_r_def, 0.932921_r_def, 0.932759_r_def, 0.932597_r_def, 0.932597_r_def, &
       0.932436_r_def, 0.932436_r_def, 0.932113_r_def, 0.932113_r_def, 0.931628_r_def, &
       0.931951_r_def, 0.931466_r_def, 0.931466_r_def, 0.931143_r_def, 0.930981_r_def, &
       0.930819_r_def, 0.930981_r_def, 0.930496_r_def, 0.930335_r_def, 0.930173_r_def, &
       0.930011_r_def, 0.930011_r_def, 0.929526_r_def, 0.929526_r_def, 0.929203_r_def, &
       0.928880_r_def, 0.928557_r_def, 0.928395_r_def, 0.928233_r_def, 0.927910_r_def, &
       0.927587_r_def, 0.927425_r_def, 0.927263_r_def, 0.926940_r_def, 0.926617_r_def, &
       0.926294_r_def, 0.926132_r_def, 0.925970_r_def, 0.925647_r_def, 0.925324_r_def, &
       0.925162_r_def, 0.924839_r_def, 0.924677_r_def, 0.924354_r_def, 0.923869_r_def, &
       0.923384_r_def, 0.923061_r_def, 0.922738_r_def, 0.922414_r_def, 0.922091_r_def, &
       0.921606_r_def, 0.921283_r_def, 0.920475_r_def, 0.920151_r_def, 0.919505_r_def, &
       0.918697_r_def /)

  ! Emissivity spectrum for sea (using the IREMIS model)
  real(r_def), parameter :: emis_sea_iremis(1:n_val) = &
    (/ 0.829415_r_def, 0.829415_r_def, 0.829415_r_def, 0.829415_r_def, 0.829415_r_def, &
       0.829415_r_def, 0.829415_r_def, 0.829415_r_def, 0.829415_r_def, 0.829415_r_def, &
       0.829415_r_def, 0.832103_r_def, 0.834120_r_def, 0.835956_r_def, 0.837573_r_def, &
       0.840842_r_def, 0.844644_r_def, 0.849024_r_def, 0.855006_r_def, 0.861022_r_def, &
       0.867069_r_def, 0.874420_r_def, 0.879676_r_def, 0.884210_r_def, 0.888080_r_def, &
       0.890979_r_def, 0.893071_r_def, 0.894495_r_def, 0.895447_r_def, 0.895840_r_def, &
       0.895818_r_def, 0.895517_r_def, 0.895083_r_def, 0.894606_r_def, 0.894128_r_def, &
       0.893755_r_def, 0.893505_r_def, 0.893316_r_def, 0.893276_r_def, 0.893305_r_def, &
       0.893386_r_def, 0.893647_r_def, 0.893866_r_def, 0.894014_r_def, 0.894171_r_def, &
       0.894333_r_def, 0.894427_r_def, 0.894456_r_def, 0.894444_r_def, 0.894380_r_def, &
       0.894305_r_def, 0.893618_r_def, 0.893433_r_def, 0.893582_r_def, 0.893988_r_def, &
       0.894563_r_def, 0.895226_r_def, 0.895995_r_def, 0.897038_r_def, 0.898089_r_def, &
       0.899133_r_def, 0.900329_r_def, 0.901405_r_def, 0.902553_r_def, 0.903789_r_def, &
       0.905040_r_def, 0.906317_r_def, 0.907594_r_def, 0.908893_r_def, 0.910179_r_def, &
       0.911600_r_def, 0.913013_r_def, 0.914720_r_def, 0.916655_r_def, 0.918588_r_def, &
       0.921045_r_def, 0.923549_r_def, 0.926072_r_def, 0.928952_r_def, 0.931907_r_def, &
       0.935064_r_def, 0.938517_r_def, 0.942251_r_def, 0.945496_r_def, 0.949562_r_def, &
       0.952738_r_def, 0.956244_r_def, 0.958788_r_def, 0.961387_r_def, 0.962776_r_def, &
       0.964222_r_def, 0.964933_r_def, 0.965058_r_def, 0.965097_r_def, 0.964472_r_def, &
       0.963873_r_def, 0.963094_r_def, 0.962288_r_def, 0.961251_r_def, 0.960211_r_def, &
       0.959665_r_def, 0.958646_r_def, 0.957691_r_def, 0.956619_r_def, 0.955974_r_def, &
       0.955248_r_def, 0.954795_r_def, 0.954262_r_def, 0.953881_r_def, 0.953517_r_def, &
       0.953182_r_def, 0.953011_r_def, 0.952887_r_def, 0.952798_r_def, 0.952312_r_def, &
       0.951653_r_def, 0.951068_r_def, 0.950612_r_def, 0.950114_r_def, 0.949646_r_def, &
       0.949128_r_def, 0.949078_r_def, 0.948993_r_def, 0.948692_r_def, 0.948627_r_def, &
       0.948562_r_def, 0.948452_r_def, 0.948345_r_def, 0.948243_r_def, 0.948217_r_def, &
       0.948200_r_def, 0.948110_r_def, 0.948057_r_def, 0.948054_r_def, 0.947998_r_def, &
       0.947915_r_def, 0.947908_r_def, 0.947910_r_def, 0.947840_r_def, 0.947771_r_def, &
       0.947609_r_def, 0.947436_r_def, 0.947199_r_def, 0.946952_r_def, 0.946646_r_def, &
       0.946346_r_def, 0.946112_r_def, 0.945878_r_def, 0.945512_r_def, 0.945138_r_def, &
       0.944758_r_def, 0.944375_r_def, 0.943984_r_def, 0.943584_r_def, 0.943152_r_def, &
       0.942553_r_def, 0.941951_r_def, 0.941215_r_def, 0.940425_r_def, 0.939629_r_def, &
       0.938797_r_def, 0.937887_r_def, 0.938665_r_def, 0.939951_r_def, 0.941149_r_def, &
       0.944554_r_def, 0.947994_r_def, 0.951021_r_def, 0.953095_r_def, 0.955037_r_def, &
       0.955844_r_def, 0.955586_r_def, 0.955232_r_def, 0.954681_r_def, 0.954023_r_def, &
       0.953350_r_def, 0.952736_r_def, 0.952175_r_def, 0.951611_r_def, 0.951109_r_def, &
       0.950692_r_def, 0.950275_r_def, 0.949876_r_def, 0.949560_r_def, 0.949245_r_def, &
       0.948930_r_def, 0.948622_r_def, 0.948317_r_def, 0.948012_r_def, 0.947737_r_def, &
       0.947527_r_def, 0.947317_r_def, 0.947106_r_def, 0.946901_r_def, 0.946698_r_def, &
       0.946496_r_def, 0.946293_r_def, 0.946164_r_def, 0.946042_r_def, 0.945920_r_def, &
       0.945798_r_def, 0.945679_r_def, 0.945559_r_def, 0.945440_r_def, 0.945321_r_def, &
       0.945239_r_def, 0.945160_r_def, 0.945081_r_def, 0.945002_r_def, 0.944972_r_def, &
       0.944966_r_def, 0.944961_r_def, 0.944955_r_def, 0.944950_r_def, 0.944948_r_def, &
       0.944947_r_def, 0.944945_r_def, 0.944943_r_def, 0.944905_r_def, 0.944844_r_def, &
       0.944782_r_def, 0.944721_r_def, 0.944659_r_def, 0.944602_r_def, 0.944546_r_def, &
       0.944489_r_def, 0.944432_r_def, 0.944375_r_def, 0.944278_r_def, 0.944167_r_def, &
       0.944055_r_def, 0.943943_r_def, 0.943831_r_def, 0.943721_r_def, 0.943614_r_def, &
       0.943507_r_def, 0.943400_r_def, 0.943293_r_def, 0.943186_r_def, 0.943083_r_def, &
       0.942980_r_def, 0.942878_r_def, 0.942775_r_def, 0.942673_r_def, 0.942568_r_def, &
       0.942446_r_def, 0.942323_r_def, 0.942200_r_def, 0.942078_r_def, 0.941955_r_def, &
       0.941832_r_def, 0.941692_r_def, 0.941552_r_def, 0.941412_r_def, 0.941272_r_def, &
       0.941132_r_def, 0.940992_r_def, 0.940844_r_def, 0.940691_r_def, 0.940537_r_def, &
       0.940384_r_def, 0.940230_r_def, 0.940077_r_def, 0.939924_r_def, 0.939724_r_def, &
       0.939516_r_def, 0.939308_r_def, 0.939100_r_def, 0.938892_r_def, 0.938684_r_def, &
       0.938477_r_def, 0.938263_r_def, 0.938048_r_def, 0.937832_r_def, 0.937617_r_def, &
       0.937401_r_def, 0.937186_r_def, 0.936971_r_def, 0.936742_r_def, 0.936464_r_def, &
       0.936187_r_def, 0.935910_r_def, 0.935632_r_def, 0.935355_r_def, 0.935078_r_def, &
       0.934800_r_def, 0.934501_r_def, 0.934145_r_def, 0.933790_r_def, 0.933434_r_def, &
       0.933081_r_def, 0.932740_r_def, 0.932397_r_def, 0.932054_r_def, 0.931711_r_def, &
       0.931319_r_def, 0.930920_r_def, 0.930521_r_def, 0.930121_r_def, 0.929628_r_def, &
       0.929038_r_def /)

  ! Planck function for tile temperature
  do i = 1, n_val
    planck_wavelentbl(i) = planck(tile_temp, wavelentbl(i))
  end do

  ! Planck weighted contribution of k'th lookup table position to j'th band
  ! and total weight for j'th band
  weight = 0.0_r_def
  total_weight = 0.0_r_def
  do j = 1, n_bands
    do k = 1, n_val
      weight(k, j) = planck_wavelentbl(k) &
            * (min(wavelen_high(j), wavelentbl_high(k)) - max(wavelen_low(j), wavelentbl_low(k)))
      if (weight(k, j) < 0.0_r_def) weight(k, j) = 0.0_r_def
    end do
    do k = 1, n_val
      total_weight(j) = total_weight(j) + weight(k, j)
    end do
  end do

  ! Correct weights for excluded bands
  do j = 1, n_bands
    do i = 1, n_band_exclude(j)
      do k = 1, n_val
        weight(k, j) = weight(k, j) - weight(k, index_band_exclude(i, j))
      end do
    total_weight(j) = total_weight(j) - total_weight(index_band_exclude(i, j))
    end do
  end do

  ! Weights for grey emissivity (each lookup table postion and total)
  grey_weight = 0.0_r_def
  grey_total_weight = 0.0_r_def
  do k = 1, n_val
    grey_weight(k) = planck_wavelentbl(k) * (wavelentbl_high(k) - wavelentbl_low(k))
    grey_total_weight = grey_total_weight + grey_weight(k)
  end do


  ! Write out band emissivities and grey emissivity for one of the data sources

  emis = 0.0_r_def
  grey_emis = 0.0_r_def

  if (surf_type == 'desert_feldman') then
    do j = 1, n_bands
      do k = 1, n_val
        emis(j) = emis(j) + weight(k, j) * emis_desert_feldman(k) / total_weight(j)
      end do
    end do
    do k = 1, n_val
      grey_emis = grey_emis + grey_weight(k) * emis_desert_feldman(k) / grey_total_weight
    end do
  end if

  if (surf_type == 'sea_feldman') then
    do j = 1, n_bands
      do k = 1, n_val
        emis(j) = emis(j) + weight(k, j) * emis_sea_feldman(k) / total_weight(j)
      end do
    end do
    do k = 1, n_val
      grey_emis = grey_emis + grey_weight(k) * emis_sea_feldman(k) / grey_total_weight
    end do
  end if

  if (surf_type == 'sea_iremis') then
    do j = 1, n_bands
      do k = 1, n_val
        emis(j) = emis(j) + weight(k, j) * emis_sea_iremis(k) / total_weight(j)
      end do
    end do
    do k = 1, n_val
      grey_emis = grey_emis + grey_weight(k) * emis_sea_iremis(k) / grey_total_weight
    end do
  end if

end subroutine specemis

end module specemis_mod
