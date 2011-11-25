(ns foobar.core
  (:require [clojure.string :as str]))

;(load-file "src/foobar/vcf.clj")

(use 'foobar.vcf)
(use '(incanter core stats charts ) )

; processing commandline for filename
(def vcf_file (nth *command-line-args* 0))


; plot QUAL
;(def a (load-vcf vcf_file))
;(with-data a
;      (view (histogram ($ :QUAL))))


; Todo
; plot AF from INFO field
; INFO-AF
;(vcf2tsv vcf_file ) this is so great!
; error can't convert String to Inte
; Integer.parseInt
; A bug with 1.000 haven't fixed yet.

(def a (load-vcf vcf_file))
(def b (with-data a ($ :INFO-AF)))

(println (str/join "\n" (map #(str/join "\t" %) b)))

;(with-data a
;     (view (histogram (map #(Double/parseDouble %) ($ :INFO-AF)))))

