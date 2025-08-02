(defun ceiling (x / i)
  ;; Просто реализиране на математическата функция ceiling
  (setq i (fix x))
  (if (< x i)
    (+ i 1)
    i))

(defun c:ContourLaplace (/ rows cols ss n p1 p2 dx dy grid fixed k ent entData pt
                            x y z row col i j tolerance maxIter iter diff
                            left right up down newVal curVal delta
                            contourInterval minZ maxZ start end levels lev
                            points z1 z2 z3 z4 xi yi)

  (setq rows 100)
  (setq cols 100)

  ;; Селектиране на точки и блокове
  (princ "\nИзберете точки (POINT) или блокове (INSERT): ")
  (setq ss (ssget '((0 . "POINT,INSERT"))))
  (if (not ss)
    (progn
      (princ "\nНе са избрани точки.")
      (princ)
    )
    (progn
      (setq n (sslength ss))

      ;; Задаване граници на мрежата
      (setq p1 (getpoint "\nПосочете долен ляв ъгъл на мрежата: "))
      (setq p2 (getpoint p1 "\nПосочете горен десен ъгъл на мрежата: "))

      ;; Изчисляване на стъпката по X и Y
      (setq dx (/ (- (car p2) (car p1)) (float (1- cols))))
      (setq dy (/ (- (cadr p2) (cadr p1)) (float (1- rows))))

      ;; Създаване на двумерни масиви (SafeArray) за височините и фиксираните клетки
      (setq grid (vlax-make-safearray vlax-vbDouble (cons 0 (1- rows)) (cons 0 (1- cols))))
      (setq fixed (vlax-make-safearray vlax-vbInteger (cons 0 (1- rows)) (cons 0 (1- cols))))

      ;; Инициализация на масивите
      (setq i 0)
      (while (< i rows)
        (setq j 0)
        (while (< j cols)
          (vlax-safearray-put-element grid (list i j) 0.0)
          (vlax-safearray-put-element fixed (list i j) 0) ; 0 = не е фиксирано
          (setq j (1+ j))
        )
        (setq i (1+ i))
      )

      ;; Задаване на височини за фиксираните точки
      (setq k 0)
      (while (< k n)
        (setq ent (ssname ss k))
        (setq entData (entget ent))
        (if (= (cdr (assoc 0 entData)) "POINT")
          (setq pt (cdr (assoc 10 entData))) ; координати на POINT
          (setq pt (cdr (assoc 10 entData))) ; базова точка на блок INSERT
        )
        (setq x (car pt))
        (setq y (cadr pt))
        (setq z (caddr pt))

        ;; Определяне на индексите в мрежата
        (setq col (fix (/ (- x (car p1)) dx)))
        (setq row (fix (/ (- y (cadr p1)) dy)))

        ;; Ограничаване в рамките на масива
        (if (< row 0) (setq row 0))
        (if (>= row rows) (setq row (1- rows)))
        (if (< col 0) (setq col 0))
        (if (>= col cols) (setq col (1- cols)))

        ;; Записване на височината и маркиране на клетката като фиксирана
        (vlax-safearray-put-element grid (list row col) z)
        (vlax-safearray-put-element fixed (list row col) 1)

        (setq k (1+ k))
      )

      ;; Решаване на дискретното уравнение на Лаплас: средно от съседните точки:contentReference[oaicite:3]{index=3}
      (setq tolerance 0.001)
      (setq maxIter 500)
      (setq iter 0)
      (setq diff tolerance)

      (while (and (< iter maxIter) (> diff tolerance))
        (setq diff 0.0)
        (setq i 0)
        (while (< i rows)
          (setq j 0)
          (while (< j cols)
            (if (= (vlax-safearray-get-element fixed (list i j)) 0)
              (progn
                ;; Изчисляване на средното от четирите съседа (със стойността на текущата клетка по границите)
                (setq left  (if (> j 0)        (vlax-safearray-get-element grid (list i (1- j))) (vlax-safearray-get-element grid (list i j))))
                (setq right (if (< j (1- cols))(vlax-safearray-get-element grid (list i (1+ j))) (vlax-safearray-get-element grid (list i j))))
                (setq up    (if (< i (1- rows))(vlax-safearray-get-element grid (list (1+ i) j)) (vlax-safearray-get-element grid (list i j))))
                (setq down  (if (> i 0)        (vlax-safearray-get-element grid (list (1- i) j)) (vlax-safearray-get-element grid (list i j))))
                (setq newVal (/ (+ left right up down) 4.0))
                (setq curVal (vlax-safearray-get-element grid (list i j)))
                (setq delta (abs (- newVal curVal)))
                (if (> delta diff) (setq diff delta))
                (vlax-safearray-put-element grid (list i j) newVal)
              )
            )
            (setq j (1+ j))
          )
          (setq i (1+ i))
        )
        (setq iter (1+ iter))
      )

      ;; Създаване на 3D мрежов обект (3DMESH)
      (command "_.3DMESH" rows cols)
      (setq i 0)
      (while (< i rows)
        (setq j 0)
        (while (< j cols)
          (setq x (+ (car p1) (* j dx)))
          (setq y (+ (cadr p1) (* i dy)))
          (setq z (vlax-safearray-get-element grid (list i j)))
          (command (list x y z))
          (setq j (1+ j))
        )
        (setq i (1+ i))
      )
      (command "") ; край на 3DMESH

      ;; Въвеждане на стъпка между контурните линии
      (setq contourInterval (getreal "\nВъведете стъпка между хоризонталите (например 1.0): "))
      (if (null contourInterval) (setq contourInterval 1.0))

      ;; Определяне на минимална и максимална височина
      (setq minZ nil)
      (setq maxZ nil)
      (setq i 0)
      (while (< i rows)
        (setq j 0)
        (while (< j cols)
          (setq z (vlax-safearray-get-element grid (list i j)))
          (if (or (null minZ) (< z minZ)) (setq minZ z))
          (if (or (null maxZ) (> z maxZ)) (setq maxZ z))
          (setq j (1+ j))
        )
        (setq i (1+ i))
      )

      ;; Формиране на списък от нива за изолиниите
      (setq levels nil)
      (setq start (fix (/ minZ contourInterval)))
      (setq end   (fix (/ maxZ contourInterval)))
      (setq i start)
      (while (<= i end)
        (setq lev (* i contourInterval))
        (setq levels (cons lev levels))
        (setq i (1+ i))
      )
      (setq levels (reverse levels))

      ;; Генериране на контурните линии чрез алгоритъма Marching Squares:contentReference[oaicite:4]{index=4}
      (foreach lev levels
        (setq i 0)
        (while (< i (1- rows))
          (setq j 0)
          (while (< j (1- cols))
            ;; Височини в четирите ъгъла на клетката
            (setq z1 (vlax-safearray-get-element grid (list i j)))
            (setq z2 (vlax-safearray-get-element grid (list i (1+ j))))
            (setq z3 (vlax-safearray-get-element grid (list (1+ i) j)))
            (setq z4 (vlax-safearray-get-element grid (list (1+ i) (1+ j))))
            (setq points nil)

            ;; проверка за пресичане по долния ръб (z1-z2)
            (if (or (and (<= z1 lev) (> z2 lev)) (and (> z1 lev) (<= z2 lev)))
              (progn
                (setq t (/ (- lev z1) (- z2 z1)))
                (setq xi (+ (+ (car p1) (* j dx)) (* t dx)))
                (setq yi (+ (cadr p1) (* i dy)))
                (setq points (cons (list xi yi lev) points))
              )
            )
            ;; десен ръб (z2-z4)
            (if (or (and (<= z2 lev) (> z4 lev)) (and (> z2 lev) (<= z4 lev)))
              (progn
                (setq t (/ (- lev z2) (- z4 z2)))
                (setq xi (+ (car p1) (* (+ j 1) dx)))
                (setq yi (+ (+ (cadr p1) (* i dy)) (* t dy)))
                (setq points (cons (list xi yi lev) points))
              )
            )
            ;; горен ръб (z3-z4)
            (if (or (and (<= z3 lev) (> z4 lev)) (and (> z3 lev) (<= z4 lev)))
              (progn
                (setq t (/ (- lev z3) (- z4 z3)))
                (setq xi (+ (+ (car p1) (* j dx)) (* t dx)))
                (setq yi (+ (cadr p1) (* (+ i 1) dy)))
                (setq points (cons (list xi yi lev) points))
              )
            )
            ;; ляв ръб (z1-z3)
            (if (or (and (<= z1 lev) (> z3 lev)) (and (> z1 lev) (<= z3 lev)))
              (progn
                (setq t (/ (- lev z1) (- z3 z1)))
                (setq xi (+ (car p1) (* j dx)))
                (setq yi (+ (+ (cadr p1) (* i dy)) (* t dy)))
                (setq points (cons (list xi yi lev) points))
              )
            )

            ;; Ако има две точки на пресичане – изчертаване на 3D полилиния
            (if (= (length points) 2)
              (progn
                (command "_.3Dpoly" (car points) (cadr points) "")
              )
            )

            (setq j (1+ j))
          )
          (setq i (1+ i))
        )
      )

      (princ "\nГотово.")
    )
  )
  (princ)
)