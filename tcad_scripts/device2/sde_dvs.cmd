;; =============================================================================
;;  Device 2: paper Fig. 5(b) short-channel HEMT (Nazir et al., IEEE TED 72(9), 2025)
;;  Step 6b-1: structure, mole fraction, doping, contacts and mesh only.
;;  Sentaurus Structure Editor, T-2022.03-SP2.  (2026-11-10)
;;  Built from GaN_SecII/sde_dvs.cmd (Sec. II, run n1, 2026-10-12); same mole-fraction
;;  windows, doping and mesh commands, changed only where Fig. 5(b) differs.
;;
;;  Tags: [P] paper (Fig. 5(b) / Sec. III text)   [R14] ref. [14] Fig. 2(a)/text
;;        [PC] M. S. Nazir, personal communication, 2026-11-10 ("use any existing
;;             Sentaurus deck, change the mole fraction to graded")
;;        [ASSUME] given in neither - revisit   [NUM] numerical only
;;
;;  Stack, bottom -> top [P, Fig. 5(b)]:
;;   substrate (stand-in) / Al0.04Ga0.96N buffer / GaN / graded AlGaN x = 0 -> 0.15 /
;;   Al0.27Ga0.73N barrier / SiN; gate 50 nm between two 0.6 um SiN spans;
;;   regrown n+ GaN under the source and drain contacts.
;;  Changes from the Sec. II deck:
;;   1. grade 0 -> 0.15 [P] with a separate 0.27 barrier [P] (Sec. II: 0 -> 0.27, tied)
;;   2. access: 0.6 um SiN each side of the gate [P, Fig. 5(b)] (Sec. II: 0.6 um total [R14])
;;   3. S/D: explicit regrown n+ GaN blocks with the metal on top [P, Fig. 5(b)]
;;      (Sec. II: metal post through the barrier plus an n+ Gaussian stand-in)
;;
;;  Coordinates: x vertical, pointing DOWN (substrate bottom x = 0, surface at
;;  negative x). Lower-case x... = bottom edge of a layer, upper-case X... = top.
;; =============================================================================

(sde:clear)
(sde:set-process-up-direction "+Z")

;; ---------------- layer thicknesses [um] ----------------
(define tSub      100e-3)  ; [ASSUME] stand-in; Fig. 5(b) says only "Substrate"
(define tBuffer   800e-3)  ; Al0.04Ga0.96N buffer [P, Fig. 5(b)]; 0.8 um [R14]
(define XBuf      0.04)    ; buffer Al fraction [P, Fig. 5(b)]
(define tChannel   40e-3)  ; GaN, thickness 400 A [R14] (not given in the paper)
(define tGraded     6e-3)  ; graded AlGaN, 6 nm [P, Sec. III]
(define tBarrier    9e-3)  ; Al0.27Ga0.73N, thickness 90 A [R14] (not given in the paper)
(define XGradeTop 0.15)    ; grade x = 0 (bottom) -> 0.15 (top) [P, Fig. 5(b)]
(define XBarMole  0.27)    ; barrier Al fraction [P, Fig. 5(b)]
(define tSiN       50e-3)  ; SiN 0.05 um [R14]
(define tGas      200e-3)  ; [NUM]

;; ---------------- lateral dimensions [um] ----------------
(define lContact 0.5)      ; [ASSUME] S/D length (not given), as in the Sec. II deck
(define Lsg 0.6)           ; SiN between each contact and the gate, 0.6 um [P, Fig. 5(b)]
(define Lg  0.05)          ; gate length 50 nm [P, Sec. III]
(define Lsd (+ Lsg Lg Lsg))                         ; source-drain spacing 1.25 um
(define DeviceLength (+ lContact Lsd lContact))     ; 2.25 um
(define yG0 (+ lContact Lsg))                       ; gate left edge,  1.10 um
(define yG1 (+ yG0 Lg))                             ; gate right edge, 1.15 um
(define yS1 lContact)                               ; inner edge of the source n+ block
(define yD0 (- DeviceLength lContact))              ; inner edge of the drain n+ block

;; ---------------- regrown n+ GaN [ASSUME: depth and doping not given] ----------------
(define dRegrow 20e-3)     ; [ASSUME] n+ GaN reaches 20 nm below the GaN top (below the 2DEG)
(define NplusConc 5e19)    ; [ASSUME] cm^-3, the peak used by the Sec. II n+ stand-in

;; ---------------- derived edges ----------------
(define xSub 0)            (define XSub (- xSub tSub))
(define xBuffer XSub)      (define XBuffer (- xBuffer tBuffer))
(define xChannel XBuffer)  (define XChannel (- xChannel tChannel))
(define xGraded XChannel)  (define XGraded (- xGraded tGraded))
(define xBarrier XGraded)  (define XBarrier (- xBarrier tBarrier))
(define xSiN XBarrier)     (define XSiN (- xSiN tSiN))
(define xGas XSiN)         (define XGas (- xGas tGas))
(define xNplus (+ XChannel dRegrow))                ; bottom edge of the n+ blocks

;; ---------------- stack ----------------
(sdegeo:create-rectangle (position xSub 0 0)     (position XSub DeviceLength 0)     "Oxide" "substrate")
(sdegeo:create-rectangle (position xBuffer 0 0)  (position XBuffer DeviceLength 0)  "AlGaN" "buffer")
(sdegeo:create-rectangle (position xChannel 0 0) (position XChannel DeviceLength 0) "GaN"   "channel")
(sdegeo:create-rectangle (position xGraded 0 0)  (position XGraded DeviceLength 0)  "AlGaN" "graded")
(sdegeo:create-rectangle (position xBarrier 0 0) (position XBarrier DeviceLength 0) "AlGaN" "barrier")

;; ---------------- regrown n+ GaN source/drain blocks [P, Fig. 5(b)] ----------------
;; Each block replaces barrier + grade + top dRegrow of GaN under its contact.
(sdegeo:set-default-boolean "ABA")      ; new bodies replace existing ones where they overlap
(sdegeo:create-rectangle (position xNplus 0 0)   (position XBarrier yS1 0)          "GaN" "nplus_source")
(sdegeo:create-rectangle (position xNplus yD0 0) (position XBarrier DeviceLength 0) "GaN" "nplus_drain")

;; ---------------- gate: on the barrier, through the SiN ----------------
(sdegeo:set-default-boolean "BAB")      ; new bodies are trimmed by existing ones
(sdegeo:create-rectangle (position XBarrier yG0 0) (position XSiN yG1 0) "Silver" "Con.Gate")

;; ---------------- SiN and gas ----------------
(sdegeo:create-rectangle (position xSiN 0 0) (position XSiN DeviceLength 0) "Nitride" "passivation")
(sde:separate-lumps)
(sdegeo:create-rectangle (position xGas 0 0) (position XGas DeviceLength 0) "Gas" "TopGas")
(sdegeo:set-default-boolean "ABA")

;; ---------------- Al mole fraction ----------------
;; Buffer: constant 0.04 [P] (region-based, cannot leak into other regions).
(sdedr:define-constant-profile "BufMole" "xMoleFraction" XBuf)
(sdedr:define-constant-profile-region "BufMolePlace" "BufMole" "buffer")
;; Grade + barrier, same NESTED windows as the Sec. II deck (placements add where
;; windows overlap): window k starts k*dx above the GaN and runs past the surface,
;; adding XGradeTop/N; the sum is a linear 0 -> 0.15 grade in the graded layer and
;; 0.15 - dMole/2 in the barrier. One more window over the barrier only adds
;; 0.27 - 0.15 + dMole/2, so the barrier is exactly 0.27 [P]. Laterally the windows
;; span only the AlGaN between the n+ blocks.
;; [NUM] As in Sec. II (checked in n1_msh.tdr), the bottom vertex row of the barrier
;; keeps the grade-top value: the 0.15 -> 0.27 step is spread over the first barrier
;; element (about 0.1 nm with this mesh).
(define nGradeSlices 60)                  ; [NUM] 0.1 nm steps
(define dxGrade (/ tGraded nGradeSlices))
(define dMole (/ XGradeTop nGradeSlices))
(define epsWin 1.3e-5)                    ; [NUM]
(do ((i 0 (+ i 1)))
    ((= i nGradeSlices))
  (let* ((xLo  (- xGraded (* i dxGrade)))
         (molX (if (= i 0) (/ dMole 2) dMole))
         (win  (string-append "GradeWin_"   (number->string i)))
         (prof (string-append "GradeMole_"  (number->string i)))
         (plc  (string-append "GradePlace_" (number->string i))))
    (sdedr:define-refeval-window win "Rectangle" (position xLo yS1 0) (position (- XBarrier epsWin) yD0 0))
    (sdedr:define-constant-profile prof "xMoleFraction" molX)
    (sdedr:define-constant-profile-placement plc prof win)))
(sdedr:define-refeval-window "BarTopWin" "Rectangle"
  (position (- xBarrier epsWin) yS1 0) (position (- XBarrier epsWin) yD0 0))
(sdedr:define-constant-profile "BarTopMole" "xMoleFraction" (+ (- XBarMole XGradeTop) (/ dMole 2)))
(sdedr:define-constant-profile-placement "BarTopPlace" "BarTopMole" "BarTopWin")

;; ---------------- doping ----------------
;; [ASSUME] unintentionally doped everywhere else (no doping given in [14] / paper)
(sdedr:define-constant-profile "ndop_const" "nSiliconActiveConcentration" 1e10)
(sdedr:define-constant-profile-region "ndop_buf" "ndop_const" "buffer")
(sdedr:define-constant-profile-region "ndop_ch"  "ndop_const" "channel")
(sdedr:define-constant-profile-region "ndop_gr"  "ndop_const" "graded")
(sdedr:define-constant-profile-region "ndop_bar" "ndop_const" "barrier")
;; [ASSUME] deep acceptors in the lower buffer, as in the Synopsys example deck
(sdedr:define-refeval-window "DeepWin" "Rectangle" (position xBuffer 0 0) (position (+ XBuffer 0.2) DeviceLength 0))
(sdedr:define-constant-profile "trap_prof" "DeepLevels" 5e16)
(sdedr:define-constant-profile-placement "trap_prof1" "trap_prof" "DeepWin" 0.07)
;; regrown n+ GaN: uniform [ASSUME]
(sdedr:define-constant-profile "nplus_prof" "NDopantActiveConcentration" NplusConc)
(sdedr:define-constant-profile-region "nplus_src" "nplus_prof" "nplus_source")
(sdedr:define-constant-profile-region "nplus_drn" "nplus_prof" "nplus_drain")

;; ---------------- contacts ----------------
(define thermal.reg.bot
  (sdegeo:create-rectangle (position (+ xSub 0.01) 0 0) (position xSub DeviceLength 0) "Aluminum" "Con.thermal"))
;; source/drain metal on top of the n+ blocks, through the SiN [P, Fig. 5(b)]
(define source
  (sdegeo:create-rectangle (position XBarrier 0 0) (position XSiN yS1 0) "Copper" "Cnt.Source"))
(define drain
  (sdegeo:create-rectangle (position XBarrier yD0 0) (position XSiN DeviceLength 0) "Copper" "Cnt.Drain"))
(define gate.reg
  (list (car (find-body-id (position (- XBarrier (/ tSiN 2)) (/ (+ yG0 yG1) 2) 0)))))

(sdegeo:define-contact-set "source")
(sdegeo:set-contact source "source" "remove")
(sdegeo:define-contact-set "drain")
(sdegeo:set-contact drain "drain" "remove")
(sdegeo:define-contact-set "gate")
(sdegeo:set-contact gate.reg "gate" "remove")
(sdegeo:define-contact-set "thermal")
(sdegeo:set-contact thermal.reg.bot "thermal" "remove")

;; ---------------- mesh [NUM] ----------------
(define xmin (position:x (car (bbox (get-body-list)))))
(define xmax (position:x (cdr (bbox (get-body-list)))))
(define ymin (position:y (car (bbox (get-body-list)))))
(define ymax (position:y (cdr (bbox (get-body-list)))))

(sdedr:define-refinement-window "whole" "Rectangle" (position xmin ymin 0) (position xmax ymax 0))
(sdedr:define-refinement-size "RefDef.all" (/ (- xmax xmin) 30) (/ (- ymax ymin) 60)
                                           (/ (- xmax xmin) 180) (/ (- ymax ymin) 300))
(sdedr:define-refinement-function "RefDef.all" "DopingConcentration" "MaxGradient" 1)
(sdedr:define-refinement-placement "meshwhole" "RefDef.all" "whole")
(sdedr:define-refinement-function "RefDef.all" "MaxLenInt" "All" "All" 10e-3 2 "DoubleSide")
(sdedr:define-refinement-function "RefDef.all" "MaxLenInt" "All" "Contact" 5e-3 1.5)
;; <= 0.5 nm at the 2DEG-forming interfaces (GaN_HEMT_FermiStat readme)
(sdedr:define-refinement-function "RefDef.all" "MaxLenInt" "channel" "graded" 0.5e-3 1.2 "DoubleSide" "UseRegionNames")
(sdedr:define-refinement-function "RefDef.all" "MaxLenInt" "graded" "barrier" 0.5e-3 1.2 "DoubleSide" "UseRegionNames")

;; barrier + grade + GaN + top 30 nm of the buffer, whole width: 1 nm vertical, 10 nm lateral
;; (the Sec. II comment said "top 30 nm of GaN"; the window below is what both decks build)
(sdedr:define-refinement-window "Ref.Win.Chan" "Rectangle"
  (position (+ xChannel 0.03) 0 0) (position XBarrier DeviceLength 0))
(sdedr:define-refinement-size "Ref.Def.Chan" 1e-3 10e-3 0.5e-3 5e-3)
(sdedr:define-refinement-placement "Ref.Place.Chan" "Ref.Def.Chan" "Ref.Win.Chan")
;; around the 50 nm gate: 2 nm lateral
(sdedr:define-refinement-window "Ref.Win.Gate" "Rectangle"
  (position (+ xChannel 0.03) (- yG0 0.05) 0) (position XBarrier (+ yG1 0.05) 0))
(sdedr:define-refinement-size "Ref.Def.Gate" 1e-3 2e-3 0.5e-3 1e-3)
(sdedr:define-refinement-placement "Ref.Place.Gate" "Ref.Def.Gate" "Ref.Win.Gate")
;; around the n+ / channel sidewalls, where the 2DEG meets the regrown GaN: 2 nm lateral
(sdedr:define-refinement-window "Ref.Win.SrcEdge" "Rectangle"
  (position (+ xNplus 0.01) (- yS1 0.02) 0) (position XBarrier (+ yS1 0.02) 0))
(sdedr:define-refinement-window "Ref.Win.DrnEdge" "Rectangle"
  (position (+ xNplus 0.01) (- yD0 0.02) 0) (position XBarrier (+ yD0 0.02) 0))
(sdedr:define-refinement-size "Ref.Def.Edge" 1e-3 2e-3 0.5e-3 1e-3)
(sdedr:define-refinement-placement "Ref.Place.SrcEdge" "Ref.Def.Edge" "Ref.Win.SrcEdge")
(sdedr:define-refinement-placement "Ref.Place.DrnEdge" "Ref.Def.Edge" "Ref.Win.DrnEdge")

(sde:build-mesh "n@node@")
