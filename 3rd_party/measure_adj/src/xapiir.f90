!auto generated from xapiir_sub.f, using ChatGPT

module xapiir_mod
   use, intrinsic :: iso_fortran_env, only: wp => real32, ip => int32
   implicit none
   private

   public :: xapiir_sub

contains

   !===========================
   ! Top-level driver (public)
   !===========================
   subroutine xapiir_sub(data, aproto, trbndw, a, iord, ftype, flo, fhi, ts, passes)
      real(wp),            intent(inout)        :: data(:)     ! filtered in-place
      character(len=2),    intent(in)           :: aproto       ! 'BU','BE','C1','C2'
      real(wp),            intent(in)           :: trbndw       ! for C1/C2
      real(wp),            intent(in)           :: a            ! for C1/C2
      integer(ip),         intent(in)           :: iord         ! <= 10
      character(len=2),    intent(in)           :: ftype        ! 'LP','HP','BP','BR'
      real(wp),            intent(in)           :: flo          ! Hz (ignored for LP)
      real(wp),            intent(in)           :: fhi          ! Hz (ignored for HP)
      real(wp),            intent(in)           :: ts           ! sampling interval (s)
      integer(ip),         intent(in)           :: passes       ! 1 => forward, 2 => zero-phase

      ! Second-order sections:  up to order 10 → at most 10 sections, each 3 coeffs
      real(wp)   :: sn(3*10), sd(3*10)
      integer(ip):: nsects
      logical    :: zp

      sn = 0.0_wp
      sd = 0.0_wp

      call design(iord, ftype(1:2), aproto(1:2), a, trbndw, flo, fhi, ts, sn, sd, nsects)
      zp = (passes /= 1_ip)
      call apply(data, zp, sn, sd, nsects)
   end subroutine xapiir_sub

   !===========================
   ! Apply IIR (SOS, optional zero-phase)
   !===========================
   subroutine apply(data, zp, sn, sd, nsects)
      real(wp),  intent(inout) :: data(:)
      logical,   intent(in)    :: zp
      real(wp),  intent(in)    :: sn(:)   ! packed SOS: [b0,b1,b2] per section
      real(wp),  intent(in)    :: sd(:)   ! packed SOS: [a0(=1),a1,a2] after bilinear
      integer(ip), intent(in)  :: nsects

      integer(ip) :: j, i, jptr
      real(wp)    :: x1, x2, y1, y2, b0, b1, b2, a1, a2, out

      if (nsects <= 0) return

      ! forward pass
      jptr = 1
      do j = 1, nsects
         x1 = 0.0_wp; x2 = 0.0_wp
         y1 = 0.0_wp; y2 = 0.0_wp
         b0 = sn(jptr)   ; b1 = sn(jptr+1) ; b2 = sn(jptr+2)
         ! sd(jptr) should be 1 after bilinear; keep general form but use a1,a2 only
         a1 = sd(jptr+1) ; a2 = sd(jptr+2)

         do i = 1, size(data)
            out = b0*data(i) + b1*x1 + b2*x2
            out = out - (a1*y1 + a2*y2)
            y2 = y1
            y1 = out
            x2 = x1
            x1 = data(i)
            data(i) = out
         end do
         jptr = jptr + 3
      end do

      if (.not. zp) return

      ! reverse pass (zero-phase)
      jptr = 1
      do j = 1, nsects
         x1 = 0.0_wp; x2 = 0.0_wp
         y1 = 0.0_wp; y2 = 0.0_wp
         b0 = sn(jptr)   ; b1 = sn(jptr+1) ; b2 = sn(jptr+2)
         a1 = sd(jptr+1) ; a2 = sd(jptr+2)

         do i = size(data), 1, -1
            out = b0*data(i) + b1*x1 + b2*x2
            out = out - (a1*y1 + a2*y2)
            y2 = y1
            y1 = out
            x2 = x1
            x1 = data(i)
            data(i) = out
         end do
         jptr = jptr + 3
      end do
   end subroutine apply

   !===========================
   ! Filter design (analog proto → mapping → bilinear)
   !===========================
   subroutine design(iord, ftype, aproto, a, trbndw, fl, fh, ts, sn, sd, nsects)
      integer(ip),       intent(in)  :: iord
      character(len=*),  intent(in)  :: ftype   ! 'LP','HP','BP','BR'
      character(len=*),  intent(in)  :: aproto  ! 'BU','BE','C1','C2'
      real(wp),          intent(in)  :: a, trbndw, fl, fh, ts
      real(wp),          intent(out) :: sn(:), sd(:)
      integer(ip),       intent(out) :: nsects

      character(len=3) :: rtype(10)
      complex(wp)      :: p(10), z(10)
      real(wp)         :: dcvalue, flw, fhw, omegar, eps, ripple

      sn = 0.0_wp; sd = 0.0_wp; nsects = 0

      select case (aproto(1:2))
      case ('BU')
         call buroots(p, rtype, dcvalue, nsects, iord)
      case ('BE')
         call beroots(p, rtype, dcvalue, nsects, iord)
      case ('C1')
         call chebparm(a, trbndw, iord, eps, ripple)
         call c1roots(p, rtype, dcvalue, nsects, iord, eps)
      case ('C2')
         omegar = 1.0_wp + trbndw
         call c2roots(p, z, rtype, dcvalue, nsects, iord, a, omegar)
      case default
         error stop 'DESIGN: unknown analog prototype'
      end select

      ! analog mapping
      select case (ftype(1:2))
      case ('BP')
         flw = warp(fl*ts/2.0_wp, 2.0_wp)
         fhw = warp(fh*ts/2.0_wp, 2.0_wp)
         call lptbp(p, z, rtype, dcvalue, nsects, flw, fhw, sn, sd)
      case ('BR')
         flw = warp(fl*ts/2.0_wp, 2.0_wp)
         fhw = warp(fh*ts/2.0_wp, 2.0_wp)
         call lptbr(p, z, rtype, dcvalue, nsects, flw, fhw, sn, sd)
      case ('LP')
         fhw = warp(fh*ts/2.0_wp, 2.0_wp)
         call lp(p, z, rtype, dcvalue, nsects, sn, sd)
         call cutoffs(sn, sd, nsects, fhw)
      case ('HP')
         flw = warp(fl*ts/2.0_wp, 2.0_wp)
         call lpthp(p, z, rtype, dcvalue, nsects, sn, sd)
         call cutoffs(sn, sd, nsects, flw)
      case default
         error stop 'DESIGN: unknown filter type'
      end select

      ! bilinear transform (in-place on SOS)
      call bilin2(sn, sd, nsects)
   end subroutine design

   !===========================
   ! Butterworth poles (normalized LP)
   !===========================
   subroutine buroots(p, rtype, dcvalue, nsects, iord)
      complex(wp),      intent(out) :: p(:)
      character(len=3), intent(out) :: rtype(:)
      real(wp),         intent(out) :: dcvalue
      integer(ip),      intent(out) :: nsects
      integer(ip),      intent(in)  :: iord

      integer(ip) :: half, k
      real(wp)    :: angle, pi

      pi = acos(-1.0_wp)
      half = iord/2
      nsects = 0

      ! odd order: add pole at -1
      if (2*half < iord) then
         p(1)     = cmplx(-1.0_wp, 0.0_wp, kind=wp)
         rtype(1) = 'SP'
         nsects   = 1
      end if

      do k = 1, half
         angle   = pi * (0.5_wp + real(2*k-1,wp) / real(2*iord,wp))
         nsects  = nsects + 1
         p(nsects)     = cmplx(cos(angle), sin(angle), kind=wp)
         rtype(nsects) = 'CP'
      end do

      dcvalue = 1.0_wp
   end subroutine buroots

   !===========================
   ! Bessel poles (normalized LP) – table-based
   !===========================
   subroutine beroots(p, rtype, dcvalue, nsects, iord)
      complex(wp),      intent(out) :: p(:)
      character(len=3), intent(out) :: rtype(:)
      real(wp),         intent(out) :: dcvalue
      integer(ip),      intent(out) :: nsects
      integer(ip),      intent(in)  :: iord

      ! (Values as in original; extend as needed)
      select case (iord)
      case (1)
         p(1)     = cmplx(-1.0_wp, 0.0_wp, kind=wp)
         rtype(1) = 'SP'
      case (2)
         p(1)     = cmplx(-1.1016013_wp,  0.6360098_wp, kind=wp)
         rtype(1) = 'CP'
      case (3)
         p(1)     = cmplx(-1.0474091_wp,  0.9992645_wp, kind=wp)
         rtype(1) = 'CP'
         p(2)     = cmplx(-1.3226758_wp,  0.0_wp,       kind=wp)
         rtype(2) = 'SP'
      case (4)
         p(1)     = cmplx(-0.9952088_wp,  1.2571058_wp, kind=wp)
         rtype(1) = 'CP'
         p(2)     = cmplx(-1.3700679_wp,  0.4102497_wp, kind=wp)
         rtype(2) = 'CP'
      case (5)
         p(1)     = cmplx(-0.9576766_wp,  1.4711244_wp, kind=wp)
         rtype(1) = 'CP'
         p(2)     = cmplx(-1.3808774_wp,  0.7179096_wp, kind=wp)
         rtype(2) = 'CP'
         p(3)     = cmplx(-1.5023160_wp,  0.0_wp,       kind=wp)
         rtype(3) = 'SP'
      case (6)
         p(1)     = cmplx(-0.9306565_wp,  1.6618633_wp, kind=wp)
         rtype(1) = 'CP'
         p(2)     = cmplx(-1.3818581_wp,  0.9714719_wp, kind=wp)
         rtype(2) = 'CP'
         p(3)     = cmplx(-1.5714904_wp,  0.3208964_wp, kind=wp)
         rtype(3) = 'CP'
      case (7)
         p(1)     = cmplx(-0.9098678_wp,  1.8364514_wp, kind=wp)
         rtype(1) = 'CP'
         p(2)     = cmplx(-1.3789032_wp,  1.1915667_wp, kind=wp)
         rtype(2) = 'CP'
         p(3)     = cmplx(-1.6120388_wp,  0.5892445_wp, kind=wp)
         rtype(3) = 'CP'
         p(4)     = cmplx(-1.6843682_wp,  0.0_wp,       kind=wp)
         rtype(4) = 'SP'
      case (8)
         p(1)     = cmplx(-0.8928710_wp,  1.9983286_wp, kind=wp)
         rtype(1) = 'CP'
         p(2)     = cmplx(-1.3738431_wp,  1.3883585_wp, kind=wp)
         rtype(2) = 'CP'
         p(3)     = cmplx(-1.6369417_wp,  0.8227968_wp, kind=wp)
         rtype(3) = 'CP'
         p(4)     = cmplx(-1.7574108_wp,  0.2728679_wp, kind=wp)
         rtype(4) = 'CP'
      case default
         error stop 'BEROOTS: order not supported in table'
      end select

      nsects  = iord - iord/2
      dcvalue = 1.0_wp
   end subroutine beroots

   !===========================
   ! Chebyshev params
   !===========================
   subroutine chebparm(a, trbndw, iord, eps, ripple)
      real(wp),   intent(in)  :: a, trbndw
      integer(ip),intent(in)  :: iord
      real(wp),   intent(out) :: eps, ripple
      real(wp) :: omegar, alpha, g

      omegar = 1.0_wp + trbndw
      alpha  = (omegar + sqrt(omegar*omegar - 1.0_wp)) ** iord
      g      = (alpha*alpha + 1.0_wp) / (2.0_wp*alpha)
      eps    = sqrt(a*a - 1.0_wp) / g
      ripple = 1.0_wp / sqrt(1.0_wp + eps*eps)
   end subroutine chebparm

   !===========================
   ! Chebyshev I poles
   !===========================
   subroutine c1roots(p, rtype, dcvalue, nsects, iord, eps)
      complex(wp),      intent(out) :: p(:)
      character(len=3), intent(out) :: rtype(:)
      real(wp),         intent(out) :: dcvalue
      integer(ip),      intent(out) :: nsects
      integer(ip),      intent(in)  :: iord
      real(wp),         intent(in)  :: eps

      integer(ip) :: half, i
      real(wp)    :: pi, gamma, s, c, angle, sigma, omega

      pi    = acos(-1.0_wp)
      half  = iord/2
      gamma = (1.0_wp + sqrt(1.0_wp + eps*eps)) / eps
      gamma = exp(log(gamma)/real(iord,wp))
      s = 0.5_wp*(gamma - 1.0_wp/gamma)
      c = 0.5_wp*(gamma + 1.0_wp/gamma)

      nsects = 0
      do i = 1, half
         rtype(i) = 'CP'
         angle    = real(2*i-1,wp) * pi / real(2*iord,wp)
         sigma    = -s * sin(angle)
         omega    =  c * cos(angle)
         p(i)     = cmplx(sigma, omega, kind=wp)
         nsects   = nsects + 1
      end do

      if (2*half < iord) then
         rtype(half+1) = 'SP'
         p(half+1)     = cmplx(-s, 0.0_wp, kind=wp)
         nsects        = nsects + 1
         dcvalue       = 1.0_wp
      else
         dcvalue       = 1.0_wp / sqrt(1.0_wp + eps*eps)
      end if
   end subroutine c1roots

   !===========================
   ! Chebyshev II poles & zeros
   !===========================
   subroutine c2roots(p, z, rtype, dcvalue, nsects, iord, a, omegar)
      complex(wp),      intent(out) :: p(:), z(:)
      character(len=3), intent(out) :: rtype(:)
      real(wp),         intent(out) :: dcvalue
      integer(ip),      intent(out) :: nsects
      integer(ip),      intent(in)  :: iord
      real(wp),         intent(in)  :: a, omegar

      integer(ip) :: half, i
      real(wp)    :: pi, gamma, s, c, angle, alpha, beta, denom, sigma, omega

      pi   = acos(-1.0_wp)
      half = iord/2
      gamma = a + sqrt(a*a - 1.0_wp)
      gamma = exp(log(gamma)/real(iord,wp))
      s = 0.5_wp*(gamma - 1.0_wp/gamma)
      c = 0.5_wp*(gamma + 1.0_wp/gamma)

      nsects = 0
      do i = 1, half
         rtype(i) = 'CPZ'
         angle = real(2*i-1,wp) * pi / real(2*iord,wp)
         alpha = -s * sin(angle)
         beta  =  c * cos(angle)
         denom = alpha*alpha + beta*beta
         sigma =  omegar * alpha / denom
         omega = -omegar * beta  / denom
         p(i)  = cmplx(sigma, omega, kind=wp)

         omega =  omegar / cos(angle)
         z(i)  = cmplx(0.0_wp, omega, kind=wp)

         nsects = nsects + 1
      end do

      if (2*half < iord) then
         rtype(half+1) = 'SP'
         p(half+1)     = cmplx(-omegar/s, 0.0_wp, kind=wp)
         nsects        = nsects + 1
      end if

      dcvalue = 1.0_wp
   end subroutine c2roots

   !===========================
   ! Prewarp frequency (for bilinear)
   !===========================
   pure real(wp) function warp(f, ts) result(w)
      real(wp), intent(in) :: f, ts
      real(wp) :: twopi, angle
      twopi = 2.0_wp*acos(-1.0_wp)
      angle = twopi*f*ts/2.0_wp
      w     = (2.0_wp*tan(angle)/ts)/twopi
   end function warp

   !===========================
   ! Build LP SOS from poles/zeros
   !===========================
   subroutine lp(p, z, rtype, dcvalue, nsects, sn, sd)
      complex(wp),      intent(in)  :: p(:), z(:)
      character(len=3), intent(in)  :: rtype(:)
      real(wp),         intent(in)  :: dcvalue
      integer(ip),      intent(in)  :: nsects
      real(wp),         intent(out) :: sn(:), sd(:)

      integer(ip) :: i, iptr
      real(wp)    :: scale

      iptr = 1
      do i = 1, nsects
         select case (rtype(i))
         case ('CPZ')
            scale         = real(p(i)*conjg(p(i)), wp) / real(z(i)*conjg(z(i)), wp)
            sn(iptr  )    = real(z(i)*conjg(z(i)), wp) * scale
            sn(iptr+1)    = -2.0_wp*real(z(i), wp)     * scale
            sn(iptr+2)    = 1.0_wp                     * scale
            sd(iptr  )    = real(p(i)*conjg(p(i)), wp)
            sd(iptr+1)    = -2.0_wp*real(p(i), wp)
            sd(iptr+2)    = 1.0_wp
         case ('CP')
            scale         = real(p(i)*conjg(p(i)), wp)
            sn(iptr  )    = scale
            sn(iptr+1)    = 0.0_wp
            sn(iptr+2)    = 0.0_wp
            sd(iptr  )    = real(p(i)*conjg(p(i)), wp)
            sd(iptr+1)    = -2.0_wp*real(p(i), wp)
            sd(iptr+2)    = 1.0_wp
         case ('SP')
            scale         = -real(p(i), wp)
            sn(iptr  )    = scale
            sn(iptr+1)    = 0.0_wp
            sn(iptr+2)    = 0.0_wp
            sd(iptr  )    = -real(p(i), wp)
            sd(iptr+1)    = 1.0_wp
            sd(iptr+2)    = 0.0_wp
         end select
         iptr = iptr + 3
      end do

      ! DC gain normalization on first section
      sn(1:3) = dcvalue * sn(1:3)
   end subroutine lp

   !===========================
   ! LP → BP analog transform (SOS out)
   !===========================
   subroutine lptbp(p, z, rtype, dcvalue, nsects, fl, fh, sn, sd)
      complex(wp),      intent(in)  :: p(:), z(:)
      character(len=3), intent(in)  :: rtype(:)
      real(wp),         intent(in)  :: dcvalue, fl, fh
      integer(ip),      intent(inout) :: nsects
      real(wp),         intent(out) :: sn(:), sd(:)

      complex(wp) :: ctmp, p1, p2, z1, z2, s, h
      real(wp)    :: pi, twopi, aa, bb, scale
      integer(ip) :: i, iptr, n

      pi    = acos(-1.0_wp)
      twopi = 2.0_wp*pi
      aa    = twopi*twopi*fl*fh
      bb    = twopi*(fh - fl)

      n = nsects
      nsects = 0
      iptr = 1

      do i = 1, n
         select case (rtype(i))
         case ('CPZ')
            ctmp = (bb*z(i))**2 - 4.0_wp*aa
            ctmp = csqrt(ctmp)
            z1   = 0.5_wp*(bb*z(i) + ctmp)
            z2   = 0.5_wp*(bb*z(i) - ctmp)

            ctmp = (bb*p(i))**2 - 4.0_wp*aa
            ctmp = csqrt(ctmp)
            p1   = 0.5_wp*(bb*p(i) + ctmp)
            p2   = 0.5_wp*(bb*p(i) - ctmp)

            ! section 1
            sn(iptr  ) = real(z1*conjg(z1), wp)
            sn(iptr+1) = -2.0_wp*real(z1, wp)
            sn(iptr+2) = 1.0_wp
            sd(iptr  ) = real(p1*conjg(p1), wp)
            sd(iptr+1) = -2.0_wp*real(p1, wp)
            sd(iptr+2) = 1.0_wp
            iptr = iptr + 3

            ! section 2
            sn(iptr  ) = real(z2*conjg(z2), wp)
            sn(iptr+1) = -2.0_wp*real(z2, wp)
            sn(iptr+2) = 1.0_wp
            sd(iptr  ) = real(p2*conjg(p2), wp)
            sd(iptr+1) = -2.0_wp*real(p2, wp)
            sd(iptr+2) = 1.0_wp
            iptr = iptr + 3

            nsects = nsects + 2

         case ('CP')
            ctmp = (bb*p(i))**2 - 4.0_wp*aa
            ctmp = csqrt(ctmp)
            p1   = 0.5_wp*(bb*p(i) + ctmp)
            p2   = 0.5_wp*(bb*p(i) - ctmp)

            sn(iptr  ) = 0.0_wp
            sn(iptr+1) = bb
            sn(iptr+2) = 0.0_wp
            sd(iptr  ) = real(p1*conjg(p1), wp)
            sd(iptr+1) = -2.0_wp*real(p1, wp)
            sd(iptr+2) = 1.0_wp
            iptr       = iptr + 3

            sn(iptr  ) = 0.0_wp
            sn(iptr+1) = bb
            sn(iptr+2) = 0.0_wp
            sd(iptr  ) = real(p2*conjg(p2), wp)
            sd(iptr+1) = -2.0_wp*real(p2, wp)
            sd(iptr+2) = 1.0_wp
            iptr       = iptr + 3

            nsects = nsects + 2

         case ('SP')
            sn(iptr  ) = 0.0_wp
            sn(iptr+1) = bb
            sn(iptr+2) = 0.0_wp
            sd(iptr  ) = aa
            sd(iptr+1) = -bb*real(p(i), wp)
            sd(iptr+2) = 1.0_wp
            iptr       = iptr + 3

            nsects = nsects + 1
         end select
      end do

      ! Gain scaling at sqrt(fl*fh)
      s = cmplx(0.0_wp, sqrt(aa), kind=wp)
      h = cmplx(1.0_wp, 0.0_wp, kind=wp)
      iptr = 1
      do i = 1, nsects
         h = h * (((sn(iptr+2)*s + sn(iptr+1))*s + sn(iptr)) / ((sd(iptr+2)*s + sd(iptr+1))*s + sd(iptr)))
         iptr = iptr + 3
      end do
      scale = dcvalue / sqrt(real(h*conjg(h), wp))
      sn(1:3) = sn(1:3) * scale
   end subroutine lptbp

   !===========================
   ! LP → BR analog transform (SOS out)
   !===========================
   subroutine lptbr(p, z, rtype, dcvalue, nsects, fl, fh, sn, sd)
      complex(wp),      intent(in)  :: p(:), z(:)
      character(len=3), intent(in)  :: rtype(:)
      real(wp),         intent(in)  :: dcvalue, fl, fh
      integer(ip),      intent(inout) :: nsects
      real(wp),         intent(out) :: sn(:), sd(:)

      complex(wp) :: cinv, ctmp, p1, p2, z1, z2
      real(wp)    :: pi, twopi, aa, bb, scale, h
      integer(ip) :: n, iptr, i

      pi    = acos(-1.0_wp)
      twopi = 2.0_wp*pi
      aa    = twopi*twopi*fl*fh
      bb    = twopi*(fh - fl)

      n = nsects
      nsects = 0
      iptr = 1

      do i = 1, n
         select case (rtype(i))
         case ('CPZ')
            cinv = 1.0_wp / z(i)
            ctmp = (bb*cinv)**2 - 4.0_wp*aa
            ctmp = csqrt(ctmp)
            z1   = 0.5_wp*(bb*cinv + ctmp)
            z2   = 0.5_wp*(bb*cinv - ctmp)

            cinv = 1.0_wp / p(i)
            ctmp = (bb*cinv)**2 - 4.0_wp*aa
            ctmp = csqrt(ctmp)
            p1   = 0.5_wp*(bb*cinv + ctmp)
            p2   = 0.5_wp*(bb*cinv - ctmp)

            sn(iptr  ) = real(z1*conjg(z1), wp)
            sn(iptr+1) = -2.0_wp*real(z1, wp)
            sn(iptr+2) = 1.0_wp
            sd(iptr  ) = real(p1*conjg(p1), wp)
            sd(iptr+1) = -2.0_wp*real(p1, wp)
            sd(iptr+2) = 1.0_wp
            iptr       = iptr + 3

            sn(iptr  ) = real(z2*conjg(z2), wp)
            sn(iptr+1) = -2.0_wp*real(z2, wp)
            sn(iptr+2) = 1.0_wp
            sd(iptr  ) = real(p2*conjg(p2), wp)
            sd(iptr+1) = -2.0_wp*real(p2, wp)
            sd(iptr+2) = 1.0_wp
            iptr       = iptr + 3

            nsects = nsects + 2

         case ('CP')
            cinv = 1.0_wp / p(i)
            ctmp = (bb*cinv)**2 - 4.0_wp*aa
            ctmp = csqrt(ctmp)
            p1   = 0.5_wp*(bb*cinv + ctmp)
            p2   = 0.5_wp*(bb*cinv - ctmp)

            sn(iptr  ) = aa
            sn(iptr+1) = 0.0_wp
            sn(iptr+2) = 1.0_wp
            sd(iptr  ) = real(p1*conjg(p1), wp)
            sd(iptr+1) = -2.0_wp*real(p1, wp)
            sd(iptr+2) = 1.0_wp
            iptr       = iptr + 3

            sn(iptr  ) = aa
            sn(iptr+1) = 0.0_wp
            sn(iptr+2) = 1.0_wp
            sd(iptr  ) = real(p2*conjg(p2), wp)
            sd(iptr+1) = -2.0_wp*real(p2, wp)
            sd(iptr+2) = 1.0_wp
            iptr       = iptr + 3

            nsects = nsects + 2

         case ('SP')
            sn(iptr  ) = aa
            sn(iptr+1) = 0.0_wp
            sn(iptr+2) = 1.0_wp
            sd(iptr  ) = -aa*real(p(i), wp)
            sd(iptr+1) =  bb
            sd(iptr+2) = -real(p(i), wp)
            iptr       = iptr + 3

            nsects = nsects + 1
         end select
      end do

      ! DC scaling: BR(dc) = LP(dc)
      h = 1.0_wp
      iptr = 1
      do i = 1, nsects
         h = h * (sn(iptr) / sd(iptr))
         iptr = iptr + 3
      end do
      scale = dcvalue / abs(h)
      sn(1:3) = sn(1:3) * scale
   end subroutine lptbr

   !===========================
   ! LP → HP analog transform (SOS out)
   !===========================
   subroutine lpthp(p, z, rtype, dcvalue, nsects, sn, sd)
      complex(wp),      intent(in)  :: p(:), z(:)
      character(len=3), intent(in)  :: rtype(:)
      real(wp),         intent(in)  :: dcvalue
      integer(ip),      intent(in)  :: nsects
      real(wp),         intent(out) :: sn(:), sd(:)

      integer(ip) :: i, iptr
      real(wp)    :: scale

      iptr = 1
      do i = 1, nsects
         select case (rtype(i))
         case ('CPZ')
            scale      = real(p(i)*conjg(p(i)), wp) / real(z(i)*conjg(z(i)), wp)
            sn(iptr  ) = 1.0_wp * scale
            sn(iptr+1) = -2.0_wp*real(z(i), wp) * scale
            sn(iptr+2) = real(z(i)*conjg(z(i)), wp) * scale
            sd(iptr  ) = 1.0_wp
            sd(iptr+1) = -2.0_wp*real(p(i), wp)
            sd(iptr+2) = real(p(i)*conjg(p(i)), wp)
         case ('CP')
            scale      = real(p(i)*conjg(p(i)), wp)
            sn(iptr  ) = 0.0_wp
            sn(iptr+1) = 0.0_wp
            sn(iptr+2) = scale
            sd(iptr  ) = 1.0_wp
            sd(iptr+1) = -2.0_wp*real(p(i), wp)
            sd(iptr+2) = real(p(i)*conjg(p(i)), wp)
         case ('SP')
            scale      = -real(p(i), wp)
            sn(iptr  ) = 0.0_wp
            sn(iptr+1) = scale
            sn(iptr+2) = 0.0_wp
            sd(iptr  ) = 1.0_wp
            sd(iptr+1) = -real(p(i), wp)
            sd(iptr+2) = 0.0_wp
         end select
         iptr = iptr + 3
      end do

      sn(1:3) = sn(1:3) * dcvalue
   end subroutine lpthp

   !===========================
   ! Adjust LP/HP cutoff (SOS)
   !===========================
   subroutine cutoffs(sn, sd, nsects, f)
      real(wp),   intent(inout) :: sn(:), sd(:)
      integer(ip),intent(in)    :: nsects
      real(wp),   intent(in)    :: f

      integer(ip) :: i, iptr
      real(wp)    :: scale

      scale = 2.0_wp*acos(-1.0_wp)*f
      iptr = 1
      do i = 1, nsects
         sn(iptr+1) = sn(iptr+1) / scale
         sn(iptr+2) = sn(iptr+2) / (scale*scale)
         sd(iptr+1) = sd(iptr+1) / scale
         sd(iptr+2) = sd(iptr+2) / (scale*scale)
         iptr = iptr + 3
      end do
   end subroutine cutoffs

   !===========================
   ! Bilinear transform (SOS in-place)
   !===========================
   subroutine bilin2(sn, sd, nsects)
      real(wp),   intent(inout) :: sn(:), sd(:)
      integer(ip),intent(in)    :: nsects

      integer(ip) :: i, iptr
      real(wp)    :: a0, a1, a2, scale

      iptr = 1
      do i = 1, nsects
         a0 = sd(iptr)
         a1 = sd(iptr+1)
         a2 = sd(iptr+2)

         scale      = a2 + a1 + a0
         sd(iptr  ) = 1.0_wp
         sd(iptr+1) = 2.0_wp*(a0 - a2) / scale
         sd(iptr+2) = (a2 - a1 + a0) / scale

         a0 = sn(iptr)
         a1 = sn(iptr+1)
         a2 = sn(iptr+2)

         sn(iptr  ) = (a2 + a1 + a0) / scale
         sn(iptr+1) = 2.0_wp*(a0 - a2) / scale
         sn(iptr+2) = (a2 - a1 + a0) / scale

         iptr = iptr + 3
      end do
   end subroutine bilin2

end module xapiir_mod


SUBROUTINE XAPIIR( DATA, NSAMPS, APROTO, TRBNDW, A, IORD, &
                  TYPE, FLO, FHI, TS, PASSES)
   
   use xapiir_mod,only : xapiir_sub
   IMPLICIT NONE
!       
   CHARACTER(len=2) ::  TYPE, APROTO                                        
   INTEGER :: NSAMPS, PASSES, IORD   
   REAL ::  DATA(NSAMPS)                                               
   REAL ::  TRBNDW, A, FLO, FHI, TS

   call xapiir_sub(data,APROTO,TRBNDW,a,iord,type,&
                flo,fhi,ts,passes)

end subroutine XAPIIR