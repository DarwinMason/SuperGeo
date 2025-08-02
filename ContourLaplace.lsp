(defun ceiling (x / i)
  ;; Просто реализиране на математическата функция ceiling
  (setq i (fix x))
  (if (< x i)
    (+ i 1)
    i)
)

(defun cross2d (a b c)
  ;; 2D cross product of vectors AB and AC
  (- (* (- (car b) (car a)) (- (cadr c) (cadr a)))
     (* (- (cadr b) (cadr a)) (- (car c) (car a))))
)

(defun convex-hull (pts / start hull point next candidate done)
  ;; Returns list of points forming the convex hull of pts (Jarvis march)
  (if (< (length pts) 3)
    pts
    (progn
      ;; find leftmost lowest point
      (setq start (car (vl-sort pts
                                '(lambda (p q)
                                   (if (= (car p) (car q))
                                     (< (cadr p) (cadr q))
                                     (< (car p) (car q)))))))
      (setq hull (list start))
      (setq point start)
      (setq done nil)
      (while (not done)
        (setq next nil)
        (foreach candidate pts
          (if (not (equal candidate point))
            (cond
              ((null next) (setq next candidate))
              ((> 0 (cross2d point next candidate)) (setq next candidate))
            )
          )
        )
        (if (or (null next) (equal next start))
          (setq done T)
          (progn
            (setq hull (cons next hull))
            (setq point next)))
      )
      (reverse hull)
    )
  )
)

(defun c:ContourLaplace (/ rows cols ss n ptList minX maxX minY maxY p1 p2 dx dy
                            grid fixed k ent entData pt x y z row col
                            i j tolerance maxIter iter diff
                            left right up down newVal curVal delta
                            contourInterval minZ maxZ start end levels lev
                            points z1 z2 z3 z4 xi yi hull)
  
  (setq rows 100
        cols 100)
  
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
      
      ;; Събиране на координати и определяне на граници
      (setq k 0 ptList nil
            minX nil maxX nil minY nil maxY nil)
      (while (< k n)
        (setq ent (ssname ss k))
        (setq entData (entget ent))
        (setq pt (cdr (assoc 10 entData)))
        (setq x (car pt) y (cadr pt) z (caddr pt))
        (setq ptList (cons pt ptList))
        (if (or (null minX) (< x minX)) (setq minX x))
        (if (or (null maxX) (> x maxX)) (setq maxX x))
        (if (or (null minY) (< y minY)) (setq minY y))
        (if (or (null maxY) (> y maxY)) (setq maxY y))
        (setq k (1+ k))
      )

      (setq p1 (list minX minY 0.0))
      (setq p2 (list maxX maxY 0.0))

      ;; Изчисляване на стъпката по X и Y
      (setq dx (/ (- maxX minX) (float (1- cols))))
      (setq dy (/ (- maxY minY) (float (1- rows))))

      ;; Създаване на двумерни масиви (SafeArray) за височините и фиксираните клетки
      (setq grid  (vlax-make-safearray vlax-vbDouble (cons 0 (1- rows)) (cons 0 (1- cols))))
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
      (foreach pt ptList
        (setq x (car pt) y (cadr pt) z (caddr pt))
        (setq col (fix (/ (- x minX) dx)))
        (setq row (fix (/ (- y minY) dy)))
        (if (< row 0) (setq row 0))
        (if (>= row rows) (setq row (1- rows)))
        (if (< col 0) (setq col 0))
        (if (>= col cols) (setq col (1- cols)))
        (vlax-safearray-put-element grid  (list row col) z)
        (vlax-safearray-put-element fixed (list row col) 1)
      )

      ;; Построяване на граница по най-външните точки
      (setq hull (convex-hull ptList))
      (if (and hull (>= (length hull) 3))
        (progn
          (command "_.PLINE")
          (foreach p hull (command p))
          (command "C" "")
        )
      )

      ;; Решаване на дискретното уравнение на Лаплас: средно от съседните точки
      (setq tolerance 0.001
            maxIter 500
            iter 0
            diff tolerance)

      (while (and (< iter maxIter) (> diff tolerance))
        (setq diff 0.0
                    i 0)
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
                      (vlax-safearray-put-element grid (list i j) newVal)))
                        (setq j (1+ j)))
                      (setq i (1+ i)))
              (setq iter (1+ iter))
      )
      
      ;; създаване на 3D мрежа
            (command "_.3DMESH" rows cols)
            (setq i 0)
            (while (< i rows)
              (setq j 0)
              (while (< j cols)
                (setq x (+ minX (* j dx))
                        y (+ minY (* i dy))
                        z (vlax-safearray-get-element grid (list i j)))
                (command (list x y z))
                (setq j (1+ j)))
                (setq i (1+ i)))
              (command "")

              ;; изчертаване на границата (конвексна обвивка)
              (setq hull (convex-hull hull))
              (if hull
                (progn
                  (command "_PLINE")
                  (foreach p hull
                    (command (append p (list 0.0))))
                  (command "C")))

              ;; изолинии
              (setq contourInterval (getreal "\nВъведете стъпка между хоризонталите (например 1.0): "))
              (if (null contourInterval) (setq contourInterval 1.0))
              (setq minZ nil
                  maxZ nil
                  i 0)
            (while (< i rows)
              (setq j 0)
              (while (< j cols)
                (setq z (vlax-safearray-get-element grid (list i j)))
                (if (or (null minZ) (< z minZ)) (setq minZ z))
                (if (or (null maxZ) (> z maxZ)) (setq maxZ z))
                (setq j (1+ j)))
              (setq i (1+ i))
            )
      
      ;; Формиране на списък от нива за изолиниите
            (setq levels nil
                  start (fix (/ minZ contourInterval))
                  end   (fix (/ maxZ contourInterval))
                  i start)
            (while (<= i end)
              (setq lev (* i contourInterval)
                    levels (cons lev levels))
              (setq i (1+ i))
            )
            (setq levels (reverse levels))

      ;; Генериране на контурните линии чрез алгоритъма Marching Squares
            (foreach lev levels
              (setq i 0)
              (while (< i (1- rows))
                (setq j 0)
                (while (< j (1- cols))
                  (setq z1 (vlax-safearray-get-element grid (list i j)))
                      (setq z2 (vlax-safearray-get-element grid (list i (1+ j))))
                      (setq z3 (vlax-safearray-get-element grid (list (1+ i) j)))
                      (setq z4 (vlax-safearray-get-element grid (list (1+ i) (1+ j))))
                      (setq points nil)
                  
                  (if (or (and (<= z1 lev) (> z2 lev)) (and (> z1 lev) (<= z2 lev)))
                    (progn
                      (setq t (/ (- lev z1) (- z2 z1)))
                      (setq xi (+ (+ minX (* j dx)) (* t dx)))
                      (setq yi (+ minY (* i dy)))
                      (setq points (cons (list xi yi lev) points))))

                  (if (or (and (<= z2 lev) (> z4 lev)) (and (> z2 lev) (<= z4 lev)))
                    (progn
                      (setq t (/ (- lev z2) (- z4 z2)))
                      (setq xi (+ minX (* (+ j 1) dx)))
                      (setq yi (+ (+ minY (* i dy)) (* t dy)))
                      (setq points (cons (list xi yi lev) points))))

                  (if (or (and (<= z3 lev) (> z4 lev)) (and (> z3 lev) (<= z4 lev)))
                    (progn
                      (setq t (/ (- lev z3) (- z4 z3)))
                      (setq xi (+ (+ minX (* j dx)) (* t dx)))
                      (setq yi (+ minY (* (+ i 1) dy)))
                      (setq points (cons (list xi yi lev) points))))

                  (if (or (and (<= z1 lev) (> z3 lev)) (and (> z1 lev) (<= z3 lev)))
                    (progn
                      (setq t (/ (- lev z1) (- z3 z1)))
                      (setq xi (+ minX (* j dx)))
                      (setq yi (+ (+ minY (* i dy)) (* t dy)))
                      (setq points (cons (list xi yi lev) points))))
                  
                  (if (= (length points) 2)
                    (command "_.3Dpoly" (car points) (cadr points) ""))

                  (setq j (1+ j)))
                (setq i (1+ i))
              )
            )
      (princ "\nГотово.")
    )
  )
)

(defun cross (o a b)
  (- (* (- (car a) (car o)) (- (cadr b) (cadr o)))
    (* (- (car b) (car o)) (- (cadr a) (cadr o)))
  )
)

(defun convex-hull (pts / hull start cur next)
  (cond
    ((< (length pts) 3) pts)
    (t
      (setq start (car (vl-sort pts
                              '(lambda (p1 p2)
                                  (if (= (car p1) (car p2))
                                      (< (cadr p1) (cadr p2))
                                      (< (car p1) (car p2)))))))
      (setq hull nil
          cur start)
      (repeat 1000
        (setq hull (append hull (list cur))
            next (car pts))
        (foreach p pts
          (if (or (= next cur) (> (cross cur next p) 0))
            (setq next p)))
        (setq cur next)
        (if (equal cur start 1e-8) (return)))
      hull))
)