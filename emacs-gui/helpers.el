(global-set-key (kbd "<f5>")
                (lambda ()
                  (interactive)
                  (save-buffer)
                  (compile (format "./libcdgbs_test %s && geo-refresh"
                                   (shell-quote-argument buffer-file-name)))))
(global-set-key (kbd "C-c h")
                (lambda (n)
                  (interactive "p")
                  (call-process "geo-load" nil 0 nil (format "h%02d.obj" n))))n
(global-set-key (kbd "C-c s")
                (lambda (n)
                  (interactive "p")
                  (call-process "geo-load" nil 0 nil (format "s%02d.obj" n))))
