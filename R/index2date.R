## translate index into date
index2date <- function(name, index, start_date, dates, ut) {
  fcase(
    name %in% c("infections", "dcases", "dab"),
    index - 1 + start_date - ut,
    name == "est_prev", dates[index],
    name == "r", index + start_date - ut,
    name == "R", index- 1 + start_date
  )
}


