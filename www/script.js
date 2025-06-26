window.addEventListener("load", function () {
  if (!localStorage.getItem("cookieConsent")) {
    const banner = document.createElement("div");
    banner.innerHTML = `
      <div id="cookie-banner" style="
        position: fixed;
        bottom: 0;
        left: 0;
        width: 100%;
        background-color: #222;
        color: white;
        text-align: center;
        padding: 15px;
        font-size: 0.95rem;
        z-index: 10000;
        box-shadow: 0 -2px 6px rgba(0,0,0,0.3);
      ">
        This site uses cookies to analyze traffic with Google Analytics.
        <a href='privacy.html' style="color: #4fd1c5; text-decoration: underline;">Learn more</a>.
        <button id="accept-cookies" style="
          margin-left: 10px;
          background-color: #4fd1c5;
          color: #000;
          border: none;
          padding: 6px 12px;
          border-radius: 3px;
          cursor: pointer;
        ">Got it</button>
      </div>
    `;
    document.body.appendChild(banner);

    document.getElementById("accept-cookies").onclick = function () {
      localStorage.setItem("cookieConsent", "true");
      document.getElementById("cookie-banner").remove();
    };
  }
});
