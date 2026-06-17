(defun highlight-stokes ()
  "highlight certain key words."
  (interactive)

  (highlight-regexp ".*linearize_top.*" 'hi-yellow)
  (highlight-regexp ".*integral_neumann_P2.*" 'hi-yellow)

  (highlight-regexp ".*linearize_bot.*" 'hi-green)
  (highlight-regexp ".*integral_robin_P2.*" 'hi-green)



  )
