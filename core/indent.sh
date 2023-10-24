
#!/bin/bash
for ext in f90 f; do
  for f in *."$ext"; do
    [ -f "$f" ] || continue
    echo "Indenting $f"
    emacs -batch "$f" --eval '(indent-region (point-min) (point-max) nil)' -f save-buffer 2> /dev/null
  done
done