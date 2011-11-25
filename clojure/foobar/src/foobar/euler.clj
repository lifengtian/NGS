
(ns foobar.euler
  (:use [foobar.utils])
  (:require [clojure.contrib.str-utils2 :only [reverse] :as str]
            [clojure.contrib.def :only [defn-memo] :as def]))

;; problem 1 - find the sum of the multiples of 3 and 5
;; less than 1000
(defn sum-of-multiples []
  (apply + (filter #(or
                     (= 0 (mod % 3))
                     (= 0 (mod % 5)))
                   (range 1 1000))))

;; problem 2 - find the sum of the even fibonacci numbers
;; less than 4,000,000
(defn fibo-seq
  ([] (concat [0 1] (fibo-seq 0 1)))
  ([a b]
     (let [n (+ a b)]
       (lazy-seq (cons n (fibo-seq b n))))))

(defn find-sum-of-even-fibo []
  (apply + (filter even?
                   (take-while #(< % 4000000) (fibo-seq)))))

;; problem 3 - find prime numbers
(def/defn-memo find-prime
  "Using the Sieve of Eratosthenes algorithm"
  [max]
  (loop [next-prime 2 result [] remaining (range 2 max)]
    (if (> (* next-prime next-prime) max)
      (concat (reverse result) remaining)
      (let [filtered (filter #(not (zero? (mod % next-prime))) remaining)]
        (recur (first filtered) (cons next-prime result) filtered)))))

(defn prime?
  "Check whether a number is prime"
  [n]
  (if (and (not (= n 2))
           (not (= n 3))
           (or (= n 1) (zero? (mod n 2)) (zero? (mod n 3))))
    false
    (every? #(not (zero? %))
            (map #(mod n %) (find-prime (inc (Math/round (Math/sqrt n))))))))

(defn truncatable-prime?
  [n]
  (loop [num n result true]
    (if (or (not result) (zero? num))
      result
      (recur (quot num 10) (prime? num)))))

(defn truncatable-prime2?
  [n]
  (loop [num n result true]
    (if (or (not result) (< num 10))
      (and (prime? num) result)
      (recur (Long/parseLong (.substring (String/valueOf num) 1)) (prime? num)))))

;; problem 4 - largest palindrom product of two three digit numbers
(defn is-palindrom? [n]
  (loop [num (String/valueOf n)]
    (cond (not (= (first num) (last num))) false
          (<= (.length num) 1) true
          :else (recur (.substring num 1 (dec (.length num)))))))

(defn find-max-palindrom-in-range [beg end]
  (reduce max
          (loop [n beg result []]
            (if (>= n end)
              result
              (recur (inc n)
                     (concat result
                             (filter #(is-palindrom? %)
                                     (map #(* n %) (range beg end)))))))))

;; problem 5 - find the smallest number divisible by the
;; numbers between 1 and 20
(defn is-divisible? [n start end]
  (every? #(= 0(mod n %)) (range start end)))

;; problem 6 - find the difference between the sum of squares
;; and the square of sums of the first 100 natural numbers
(defn sum-of-squares [n]
  (apply + (map #(* % %) (range 1 (inc n)))))

(defn square-of-sums [n]
  (let [sum (apply + (range 1 (inc n)))]
    (* sum sum)))

;; problem 9 - pythagorian triples
(defn generate-triple [n]
  (let [a (inc (* 2 n)) b (* (* 2 n) (inc n)) c (inc b)]
    [a b c (+ a b c)]))

(defn generate-triple2 [n]
  (let [a (* 2 n) b (dec (* n n)) c (inc (* n n))]
    [a b c (+ a b c)]))

(defn generate-triple3 [n]
  (loop [m (inc n)]
    (let [a (- (* m m) (* n n))
          b (* 2 (* m n)) c (+ (* m m) (* n n)) sum (+ a b c)]
      (println "a = " a "b = " b "c = " c "sum =" sum)
      (if (>= sum 1000)
        [a b c sum]
        (recur (inc m))))))

;; problem 12 - triangle number with 500 divisors
(defn find-proper-divisors [n]
  (filter #(= 0 (rem n %)) (range 1 (inc (quot n 2)))))

(defn find-divisors [n]
  (concat (find-proper-divisors n) [n]))

;; problem 14 - longest sequence
(defn gen-sequence [n]
  (loop [x n result []]
    (cond (= x 1) (reverse (cons x result))
          (even? x) (recur (/ x 2) (cons x result))
          (odd? x) (recur (inc (* x 3)) (cons x result)))))

;; problem 20 - sum of digits in 100!
(def/defn-memo factorial [n]
  (loop [n n result 1]
    (if (or (= n 0) (= n 1))
      result
      (recur (dec n) (* result n)))))

(defn sum-of-digits [n]
  (apply + (num-to-digits n)))

;; problem 21 - amicable numbers
(defn amicable? [n]
  (let [sum-of-proper-divisors (apply + (find-proper-divisors n))
        other-sum (apply + (find-proper-divisors sum-of-proper-divisors))]
    (and (not (= sum-of-proper-divisors other-sum)) (= n other-sum))))

(defn sum-of-amicable-pairs-in-range [beg end]
   (apply + (filter amicable? (range beg end))))

;; problem 22 - values of sorted names
(defn problem-22-solution []
  (apply + (map #(* (first %) (second %)) (map vector (iterate inc 1) (map word-value (map #(.substring % 1 (dec (count %))) (sort (vec (.split (slurp "/home/bozhidar/projects/clojure/project-euler/resources/names.txt") ",")))))))))

;; problem 25 - find first fibo with 1000 digits
(defn find-first-fibo-with-n-digits [n]
  (count (take-while #(< (.length (String/valueOf %)) n) (fibo-seq))))

;; problem 34 - sum of factorial
(defn sum-of-factorial [n]
  (apply + (map #(factorial (Integer/parseInt (.toString %)))
                (String/valueOf n))))

;; problem 35 - find circular primes
(defn circular? [n]
  (let [num (String/valueOf n) other (reverse num)]
    (prime? other)))

(defn prime? [n]
  (= n (last (find-prime (inc n)))))

;; problem 36
(defn dec-to-bin [n]
        (loop [num n res [] ]
           (if (= 0 num)
               (apply str res)
               (recur (quot num 2) (cons (rem num 2) res)))))

(defn find-pali-dec-and-binary [beg end]
  (apply +
         (filter #(and (is-palindrom? %) (is-palindrom? (dec-to-bin %)))
                 (range beg end))))

;; problem 41 - pandigital primes
(defn pandigital? [n]
  (.startsWith "123456789" (apply str (sort (String/valueOf n)))))

;; problem 42 - triangle words in English
(defn problem-42-solution []
  (count (filter triangle?
                 (map word-value
                      (vec (.split (slurp "/home/bozhidar/projects/clojure/project-euler/resources/words.txt") ","))))))

;; problem 45
(defn problem-45-solution []
  (last (filter
         #(and (pentagonal? %) (hexagonal? %))
         (take 100000 (triangle-numbers)))))

;; problem 48 - powers