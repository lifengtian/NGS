(ns foobar.utils
  (:require [clojure.contrib.str-utils2 :only [reverse] :as str]
            [clojure.contrib.def :only [defn-memo] :as def]))

(defn pow
  "Raise n to the power p"
  [n p]
  (apply * (repeat p n)))

(defn num-to-digits
  "Converts a number to a sequence of its digits"
  [n]
  (map #(Integer/parseInt (.toString %)) (String/valueOf n)))

(defn sum-of-digits-on-power
  "Calculates the sum of digits on a number on some power"
  [n p]
  (apply + (map #(pow % p) (num-to-digits n))))

(defn word-value
  "Calculates a words numeric value"
  [w]
  (apply + (map #(if (Character/isLetter %)
                   (- (.codePointAt (.toLowerCase (String/valueOf %)) 0) 96)
                   0) w)))

(def/defn-memo pentagonal-numbers
  ([] (pentagonal-numbers 1))
  ([n] (lazy-seq (cons (/ (* n (- (* n 3) 1)) 2)
                       (pentagonal-numbers (inc n))))))

(def/defn-memo hexagonal-numbers
  ([] (hexagonal-numbers 1))
  ([n] (lazy-seq (cons (* n (- (* n 2) 1))
                       (hexagonal-numbers (inc n))))))

(def/defn-memo triangle-numbers
  ([] (triangle-numbers 1))
  ([n] (lazy-seq (cons (/ (* n (+ n 1)) 2) (triangle-numbers (inc n))))))

(defn hexagonal? [n]
  (= n (first (drop-while #(> n %) (hexagonal-numbers)))))

(defn pentagonal? [n]
  (= n (first (drop-while #(> n %) (pentagonal-numbers)))))

(defn triangle?
  [n]
  (= n (first (drop-while #(> n %) (triangle-numbers)))))



