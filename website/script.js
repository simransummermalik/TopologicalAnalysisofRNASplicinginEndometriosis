/* Small conveniences only. All teaching content is written in the HTML. */
"use strict";

// Copy never runs a command. If the clipboard is unavailable, select the
// text so it can still be copied by hand.
document.querySelectorAll("[data-copy]").forEach(function (button) {
  button.hidden = false;
  button.addEventListener("click", async function () {
    const code = document.getElementById(button.dataset.copy);
    const status = document.getElementById("copy-status");
    if (!code || !status) return;
    try {
      if (!navigator.clipboard) throw new Error("Clipboard unavailable");
      await navigator.clipboard.writeText(code.textContent);
      status.textContent = "Copied. Paste the command into your terminal.";
    } catch (error) {
      const range = document.createRange();
      range.selectNodeContents(code);
      const selection = window.getSelection();
      selection.removeAllRanges();
      selection.addRange(range);
      status.textContent = "Command selected. Press Ctrl+C or Command+C to copy it.";
    }
  });
});

// Without JavaScript, ordinary anchors work and all steps remain visible.
const panels = Array.from(document.querySelectorAll("[data-step]"));
if (panels.length) {
  const links = Array.from(document.querySelectorAll("[data-step-link]"));
  const previous = document.getElementById("previous-step");
  const next = document.getElementById("next-step");
  const progress = document.getElementById("step-progress");
  let current = 0;

  function showStep(id, focusHeading) {
    const found = panels.findIndex(function (panel) { return panel.id === id; });
    current = found === -1 ? 0 : found;
    panels.forEach(function (panel, index) { panel.hidden = index !== current; });
    links.forEach(function (link, index) {
      if (index === current) link.setAttribute("aria-current", "step");
      else link.removeAttribute("aria-current");
    });
    previous.disabled = current === 0;
    next.disabled = current === panels.length - 1;
    progress.textContent = "Step " + (current + 1) + " of " + panels.length;
    if (focusHeading) panels[current].querySelector("h2").focus({ preventScroll: true });
  }

  function goTo(index) {
    const id = panels[index].id;
    history.pushState(null, "", "#" + id);
    showStep(id, true);
  }

  links.forEach(function (link, index) {
    link.addEventListener("click", function (event) {
      event.preventDefault();
      goTo(index);
    });
  });
  previous.addEventListener("click", function () { if (current > 0) goTo(current - 1); });
  next.addEventListener("click", function () { if (current < panels.length - 1) goTo(current + 1); });
  window.addEventListener("popstate", function () { showStep(location.hash.slice(1), false); });
  window.addEventListener("hashchange", function () { showStep(location.hash.slice(1), false); });
  document.getElementById("step-controls").hidden = false;
  showStep(location.hash.slice(1), false);
}
