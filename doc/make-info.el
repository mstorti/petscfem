(find-file "options3.texi")
(erase-buffer)
(insert-file "options.texi")
(search-forward "__INSERT_HERE__")
(beginning-of-line)
(insert-file "options2.texi")
(save-buffer)
(texinfo-all-menus-update)
(texinfo-every-node-update)
(texinfo-all-menus-update)
(save-buffer)
