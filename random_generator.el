(defun generate-random-numbers (n limit)
  (let ((seq ())
        (text ""))
    (while (< (length seq) n)
      (setq seq (cons (% (random) limit) seq)))
    (dolist (number seq text)
      (setq text
            (concat
             (number-to-string number)
             (if (> (length text) 0) ", " "") 
             text)))
    text))
    
  
(dotimes (n 30 v) (print (generate-random-numbers 10 200)))
