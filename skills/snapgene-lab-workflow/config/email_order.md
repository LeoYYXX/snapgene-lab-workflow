# Primer order email (draft only)

**Sending policy:** fill the order sheet, draft the message, and show both to the user. Do not send anything to a vendor until the user explicitly confirms. A draft is not a sent order.

Recommended attachment: a filled copy of the vendor sheet in `orders/drafts/`, never a modified template in `orders/templates/`.

## Confirmation round (to the user)

- To: `researcher@example.org` (replace with the user's address from a private local config)
- Subject: `{date}-{group}-{purpose}` (proposed vendor subject)
- Body: this is a proposed vendor email plus the filled sheet. Please reply "ok" / "确认" if I should send it to the vendor. I will not send it until you do.

## Vendor round (only after explicit confirmation)

- To: `orders@example-vendor.com` (example address; use the lab's real vendor address from private config)
- Same subject and body and the same attachment.

Do not invent a vendor contact. If the lab has not given one, stop and ask.
